## EvalState.jl
# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - EvalState
# Function to evaluate simulations of the same state
# ---
# created by Sebastian Schmitt, 30.04.2020
# ------------------------------------------------------------------------------

# Main Function
function eval_state(subfolder::Array{String,1}, inpar::Opts)

    # Get excluded folders
    what_include = Bool[]
    for subf in subfolder
        pos = findlast(isequal('/'),subf)
        if !isnothing(pos) && occursin("_EXCLUDE",subf[pos+1:end])
            append!(what_include,false)
        else
            append!(what_include,true)
        end
    end
    subfolder = subfolder[what_include]

    pos = findlast(isequal('/'),subfolder[1])-1
    folder = subfolder[1][1:pos]

    # Loading info file
    moltype, dt, natoms, molmass = load_info(subfolder[1])

    # T, p, ρ
    T, p, ρ, x, c = StaticProperties(subfolder)
    println(string("ave_state DONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS")))

    # Calculation of box length L_box
    L_box = get_L_box(subfolder[1])

    setTDM = OptsTDM(folder,subfolder,false,NaN,NaN,NaN,"","","",inpar.n_boot)
    setTDM.tskip = 2
    setTDM.cutcrit = inpar.cutcrit

    # Transport Properties
    η, η_V, D, λ = transport_properties(T.val, L_box, setTDM)
    println(string("transport_properties DONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS")))

    # Output Data
    output_results(ResultsDat(T,p,ρ,x,[],[],[],c,η,η_V,D,λ,[]),folder)
end

# Subfunctions
# Average of T, p, ρ
function StaticProperties(subf)
    # Make array
    n = length(subf)
    Tmat = zeros(n)
    pmat = zeros(n)
    ρmat = zeros(n)
    x = Array{SingleDat,1}(undef,0)
    cmat = zeros(n)
    for i = 1:n
        res = load_result(string(subf[i],"/result.dat"))
        Tmat[i] = res.T.val
        pmat[i] = res.p.val
        ρmat[i] = res.ρ.val
        if (i == 1) x = res.x end
        if typeof(res.c) == SingleDat
            cmat[i] = res.c.val
        end
    end

    # Save in structure
    T = SingleDat(mean(Tmat), std(Tmat), std(Tmat)/sqrt(n))
    p = SingleDat(mean(pmat), std(pmat), std(pmat)/sqrt(n))
    ρ = SingleDat(mean(ρmat), std(ρmat), std(ρmat)/sqrt(n))
    c = SingleDat(mean(cmat), std(cmat), std(cmat)/sqrt(n))

    return T, p, ρ, x, c
end

# Load result file
function load_result(file)
    fID = open(file,"r")
    lines = readlines(fID);
    close(fID)
    res = ResultsDat([],[],[],SingleDat[],[],[],[],[],[],[],SingleDat[],[],[])

    for i = 1:length(lines)
        pos_colon = findfirst(isequal(':'),lines[i])
        pos_open = findfirst(isequal('('),lines[i])
        pos_comma = findfirst(isequal(','),lines[i])
        pos_close = findfirst(isequal(')'),lines[i])
        if length(lines[i])>1 && !(lines[i][1] == '#')
            if !occursin("---",lines[i])
                num = 0
                if lines[i][1] == 'ρ' || lines[i][1] == 'η' || lines[i][1] == 'λ'
                    name = lines[i][1:pos_colon-2]
                else
                    name = lines[i][1:pos_colon-1]
                    if lines[i][pos_colon-1] in "123456789"
                        num = parse(Int64, lines[i][pos_colon-1])
                    end
                end
                if isnothing(pos_open) && isnothing(pos_close)
                    val = parse(Float64,strip(lines[i][pos_colon+1:pos_comma-1]))
                    std = NaN; err = NaN
                elseif pos_close < pos_comma
                    val = parse(Float64,strip(lines[i][pos_colon+1:pos_open-1]))
                    std = parse(Float64,strip(lines[i][pos_open+1:pos_close-1]))
                    err = NaN
                elseif (pos_open < pos_comma) && (pos_close > pos_comma)
                    val = parse(Float64,strip(lines[i][pos_colon+1:pos_open-1]))
                    std = parse(Float64,strip(lines[i][pos_open+1:pos_comma-1]))
                    err = parse(Float64,strip(lines[i][pos_comma+1:pos_close-1]))
                end
                if (name == "T")    res.T = SingleDat(val,std,err) end
                if (name == "p")    res.p = SingleDat(val,std,err) end
                if (name == "ρ")    res.ρ = SingleDat(val,std,err) end
                if (name == string("x",num))
                                    res.x = vcat(res.x, SingleDat(val,std,err))
                                    if length(res.x) != num @warn("res.x doesn't fit num!") end
                end
                if (name == "Etot") res.Etot = SingleDat(val,std,err) end
                if (name == "Ekin") res.Ekin = SingleDat(val,std,err) end
                if (name == "Epot") res.Epot = SingleDat(val,std,err) end
                if (name == "c")    res.c = SingleDat(val,std,err) end
                if (name == "η")    res.η = SingleDat(val,std,err) end
                if (name == "η_")   res.η_V = SingleDat(val,std,err) end
                if (name == string("D",num))
                                    res.D = vcat(res.D, SingleDat(val,std,err))
                                    if length(res.D) != num @warn("res.D doesn't fit num!") end
                end
                if (name == "λ")    res.λ = SingleDat(val,std,err) end
            end
        end
    end
    return res
end

## Function to get box length
function get_L_box(folder)
    # Get data Files
    files = readdir(folder)
    datafiles = files[startswith.(files,"data.")]

    if sum(contains.(datafiles,"NVT")) > 0
        datafile = datafiles[contains.(datafiles,"NVT")][end]
    else
        datafile = datafiles[end]
    end

    lines = readlines("$folder/$datafile")
    k = findfirst(contains.(lines,"xlo xhi"))
    L_box = parse(Float64,split(lines[k])[2]) - parse(Float64,split(lines[k])[1])

    return L_box*1e-10
end
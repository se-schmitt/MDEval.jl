# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - EvalState
# Function to evaluate simulations of the same state
# ---
# created by Sebastian Schmitt, 30.04.2020
# ------------------------------------------------------------------------------

# Main Function
function EvalState(subfolder::Array{String,1})

    inpar = read_input()

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

    # Initialization of info structure
    moltype, dt, natoms, molmass = load_info(subfolder[1])
    mass = natoms*molmass/NA    # [mass] = g

    # T, p, ρ
    T, p, ρ = StaticProperties(subfolder)
    println(string("ave_state DONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS")))

    state = state_info(T.val,p.val,ρ.val,natoms,mass)

    # Transport Properties
    η, η_V, D, λ = TransportProperties(state,set_TDM(folder,subfolder,false,2.0,0.4,NaN,"","","",inpar.n_boot))
    println(string("TransportProperties DONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS")))

    # Output Data
    OutputResult(results_struct(T,p,ρ,[],[],[],η,η_V,D,λ),folder)
end

# Subfunctions
# Average of T, p, ρ
function StaticProperties(subf)
    # Make array
    n = length(subf)
    Tmat = zeros(n)
    pmat = zeros(n)
    ρmat = zeros(n)
    for i = 1:n
        res = load_result(string(subf[i],"/result.dat"))
        Tmat[i] = res.T.val
        pmat[i] = res.p.val
        ρmat[i] = res.ρ.val
    end

    # Save in structure
    T = single_dat(mean(Tmat), std(Tmat), std(Tmat)/n)
    p = single_dat(mean(pmat), std(pmat), std(pmat)/n)
    ρ = single_dat(mean(ρmat), std(ρmat), std(ρmat)/n)

    return T, p, ρ
end

# Load result file
function load_result(file)
    fID = open(file,"r")
    lines = readlines(fID);
    close(fID)
    res = results_struct([],[],[],[],[],[],[],[],[],[])

    for i = 1:length(lines)
        pos_colon = findfirst(isequal(':'),lines[i])
        pos_open = findfirst(isequal('('),lines[i])
        pos_comma = findfirst(isequal(','),lines[i])
        pos_close = findfirst(isequal(')'),lines[i])
        if length(lines[i])>1 && !(lines[i][1] == '#')
            if !occursin("---",lines[i])
                if lines[i][1] == 'ρ' || lines[i][1] == 'η' || lines[i][1] == 'λ'
                    name = lines[i][1:pos_colon-2]
                else
                    name = lines[i][1:pos_colon-1]
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
                if (name == "T") res.T = single_dat(val,std,err) end
                if (name == "p") res.p = single_dat(val,std,err) end
                if (name == "ρ") res.ρ = single_dat(val,std,err) end
                if (name == "Etot") res.Etot = single_dat(val,std,err) end
                if (name == "Ekin") res.Ekin = single_dat(val,std,err) end
                if (name == "Epot") res.Epot = single_dat(val,std,err) end
                if (name == "η") res.η = single_dat(val,std,err) end
                if (name == "η_") res.η_V = single_dat(val,std,err) end
                if (name == "D") res.D = single_dat(val,std,err) end
                if (name == "λ") res.λ = single_dat(val,std,err) end
            end
        end
    end
    return res
end

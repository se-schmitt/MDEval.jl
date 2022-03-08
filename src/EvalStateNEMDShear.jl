## EvalStateNEMD.jl
# ------------------------------------------------------------------------------
# Evaluation Software for NEMD simulations of same state - EvalStateNEMD
# Functions to calculate Newtonian viscosity for state
# ---
# created by Jens Wagner, 29.10.2021
# ------------------------------------------------------------------------------

# Main Function
function EvalStateNEMD(subfolder::Array{String,1}, inpar::input_struct)

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

    # Get static properties
    T, p, ρ = StaticProperties(subfolder)

    # Get NEMD results of single simulations
    η_all, s_rate_all = collectNEMD(subfolder)

    # Fit η_all
    η_val = valArray(η_all)
    s_rate_val = valArray(s_rate_all)
    s_rate_Carreau, η_Carreau, p_Carreau = FitCarreau(η_val,s_rate_val)
    s_rate_Eyring, η_Eyring, p_Eyring = FitEyring(η_val,s_rate_val)

    # s_rate-η plot
    figure()
    if !(reduced_units)
        ttxt = string("Viscosity ", L"η", " at ", L"T", " = ", round(T.val, digits=2), " K, ", L"p", " = ", round(p.val, digits=1), " MPa, ", L"ρ", " = ", round(ρ.val,digits=5), " g/ml")
    elseif reduced_units
        ttxt = string("Viscosity ", L"η", " at ", L"T", " = ", round(T.val, digits=2), ", ", L"p", " = ", round(p.val, digits=1), ", ", L"ρ", " = ", round(ρ.val,digits=5))
    end
    title(ttxt)

    # Single simulation data
    loglog(s_rate_all[1].val, η_all[1].val, c = :blue, marker = "o", label = "Single simulaion")
    errorbar(s_rate_all[1].val, η_all[1].val, η_all[1].std, c = :black, elinewidth = 0.5, capsize = 1)
    for i = 2:length(η_all)
        loglog(s_rate_all[i].val, η_all[i].val, c = :blue, marker = "o")
        errorbar(s_rate_all[i].val, η_all[i].val, η_all[i].std, c = :black, elinewidth = 0.5, capsize = 1)
    end

    # Fit
    plot(s_rate_Carreau, η_Carreau, c = :orange, linestyle = "solid", label = "Carreau")
    plot(s_rate_Eyring, η_Eyring, c =:green, linestyle = "dashdot", label = "Eyring")

    if !(reduced_units)
        xlabel(string(L"γ"," / 1/s"))
        ylabel(string(L"η"," / Pa·s"))
    elseif reduced_units
        xlabel(L"γ")
        ylabel(L"η")
    end
    tick_params(axis="both",which="both",direction="in")
    legend()
    tight_layout()

    savefig(string(folder,"/viscosity.pdf"))
    close()
    Fit_file(folder, p_Carreau, p_Eyring)
end


# Subfunctions
# Collect NEMD properties from single simulations
function collectNEMD(subf)
    n = length(subf)
    η_all = Array{single_dat,1}(undef,n)
    s_rate_all = Array{single_dat,1}(undef,n)

    for i = 1:n
        res = load_result_NEMD(string(subf[i],"/result.dat"))
        η_all[i] = res.η
        s_rate_all[i] = res.s_rate
    end

    return η_all, s_rate_all
end


# Fit Carreau model
function FitCarreau(η_val, s_rate_val)
    # Function to fit η_all (p[1]: η_N, p[2]: s_rate_0, p[3]: n)
    @. fun_Carreau(x,p) = p[1]*(1+(x/p[2])^2)^((p[3]-1)/2)

    # Sort s_rate and belonging η
    η_in = copy(η_val)
    s_rate_in = copy(s_rate_val)
    s_rate_sort, η_sort = valSort(s_rate_in, η_in)

    # Fit η_all
    # Find p[2]_0 (p[2]: s_rate_0)
    idx1 = 1
    for j = 2:length(η_sort)-2
        if η_sort[j] <= 0.95 * η_sort[1] && η_sort[j+1] <= 0.95 * η_sort[1] && η_sort[j+2] <= 0.95 * η_sort[1]
            idx1 = j
            break
        end
    end

    # Create Fits for different p[3]_0 (p[3]: n)
    p3_0_Carreau = [0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5,0.45]
    fit_η = Array{Vector{Float64}}(undef, length(p3_0_Carreau))
    η_Carreau = Array{Vector{Float64}}(undef, length(p3_0_Carreau))
    s_rate_Carreau = LinRange(minimum(s_rate_val), maximum(s_rate_val), 1000)
    for j = 1:length(p3_0_Carreau)
        fit_η[j] = curve_fit(fun_Carreau, s_rate_val, η_val, [η_sort[1], s_rate_sort[idx1], p3_0_Carreau[j]]).param
        η_Carreau[j] = fun_Carreau(s_rate_Carreau, fit_η[j])
    end

    # Calculate Error of Fits
    err_η = zeros(length(p3_0_Carreau))
    for j = 1:length(p3_0_Carreau)
        for k = idx1:length(s_rate_sort)
            idx2 = findmin(abs.(s_rate_Carreau.-s_rate_sort[k]))[2]
            err_η[j] = err_η[j] + abs(η_sort[k]-η_Carreau[j][idx2])
        end
    end

    # Choose Fit with lowest Error
    idx3 = findmin(err_η)[2]
    η_N_Carreau = fit_η[idx3][1]
    s_rate_0 = fit_η[idx3][2]
    n = fit_η[idx3][3]
    η_Carreau = fun_Carreau(s_rate_Carreau, fit_η[idx3])

    return s_rate_Carreau, η_Carreau, [η_N_Carreau, s_rate_0, n]
end


# Fit Eyring model
function FitEyring(η_val,s_rate_val)
    # Function to fit η_all (p[1]: η_N, p[2]: σ_E)
    @. fun_Eyring(x,p) = p[2]/x*asinh(p[1]*x/p[2])

    # Sort s_rate and belonging η
    η_in = copy(η_val)
    s_rate_in = copy(s_rate_val)
    s_rate_sort, η_sort = valSort(s_rate_in, η_in)

    # Find s_rate_0
    idx1 = 1
    for j = 2:length(η_sort)-2
        if η_sort[j] <= 0.95 * η_sort[1] && η_sort[j+1] <= 0.95 * η_sort[1] && η_sort[j+2] <= 0.95 * η_sort[1]
            idx1 = j
            break
        end
    end

    # Create Fits for different p[2]_0 (p[2]: σ_E)
    p2_0_Eyring = [5e9,1e9,7.5e8,5e8,2.5e8,1e8,5e7,1e7]
    fit_η = Array{Vector{Float64}}(undef, length(p2_0_Eyring))
    η_Eyring = Array{Vector{Float64}}(undef, length(p2_0_Eyring))
    s_rate_Eyring = LinRange(minimum(s_rate_val), maximum(s_rate_val), 1000)
    for j = 1:length(p2_0_Eyring)
        fit_η[j] = curve_fit(fun_Eyring, s_rate_val, η_val, [η_sort[1], p2_0_Eyring[j]]).param
        η_Eyring[j] = fun_Eyring(s_rate_Eyring, fit_η[j])
    end

    # Calculate Error of Fits
    err_η = zeros(length(p2_0_Eyring))
    for j = 1:length(p2_0_Eyring)
        for k = idx1:length(s_rate_sort)
            idx2 = findmin(abs.(s_rate_Eyring.-s_rate_sort[k]))[2]
            err_η[j] = err_η[j] + abs(η_sort[k]-η_Eyring[j][idx2])
        end
    end

    # Choose Fit with lowest Error
    idx3 = findmin(err_η)[2]
    η_N_Eyring = fit_η[idx3][1]
    σ_E = fit_η[idx3][2]
    η_Eyring = fun_Eyring(s_rate_Eyring, fit_η[idx3])

    return s_rate_Eyring, η_Eyring, [η_N_Eyring, σ_E]
end


# single_dat.val to valArray
function valArray(x)
    x_val = []
    for i = 1:length(x)
        append!(x_val,x[i].val)
    end
    return x_val
end


# Sort valArray
function valSort(x_in,y_in)

    x_sort = []
    y_sort = []
    for i = 1:length(x_in)
        idx = findmin(x_in)[2]
        append!(x_sort,x_in[idx])
        append!(y_sort,y_in[idx])
        deleteat!(x_in,idx)
        deleteat!(y_in,idx)
    end
    return x_sort, y_sort
end

# Create Fit file
function Fit_file(folder, p_Carreau, p_Eyring)

    # Create Filepath
    path = string(folder,"/viscosity_fit.dat")

    # Write to file
    fID = open(path,"w")
    line1 = string("# Created by MD - Bulk Evaluation, Folder: ", folder)
    line2 = "# Format: val; [T]=K, [p]=MPa, [ρ]=g/ml, [η]=Pa*s, [s_rate]=1/s"
    if (reduced_units) line2 = "# Format: val; reduced units" end
    header = string(line1,"\n",line2,"\n")
    print(fID,header)

    # Write data
    line3 = "\n# Carreau model\n"
    print(fID,line3)
    print_prop(fID, p_Carreau[1], "η_N")
    print_prop(fID, p_Carreau[2], "γ_0")
    print_prop(fID, p_Carreau[3], "n")
    line4 = "\n# Eyring model\n"
    print(fID,line4)
    print_prop(fID, p_Eyring[1], "η_N")
    print_prop(fID, p_Eyring[2], "σ_E")
    close(fID)
end

# Load NEMD result file
function load_result_NEMD(file)
    fID = open(file,"r")
    lines = readlines(fID);
    close(fID)
    res = results_struct_nemd([],[],[],[],[],[],[],[],[],[])

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
                if (name == "T")        res.T = single_dat(val,std,err) end
                if (name == "p")        res.p = single_dat(val,std,err) end
                if (name == "ρ")        res.ρ = single_dat(val,std,err) end
                if (name == string("x",num))
                                        res.x = vcat(res.x, single_dat(val,std,err))
                                        if length(res.x) != num @warn("res.x doesn't fit num!") end
                end
                if (name == "Etot")     res.Etot = single_dat(val,std,err) end
                if (name == "Ekin")     res.Ekin = single_dat(val,std,err) end
                if (name == "Epot")     res.Epot = single_dat(val,std,err) end
                if (name == "η")        res.η = single_dat(val,std,err) end
                if (name == "s_rate")   res.s_rate = single_dat(val,std,err) end
                if (name == "λ")        res.λ = single_dat(val,std,err) end
            end
        end
    end
    return res
end

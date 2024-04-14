## EvalStateNEMD.jl
# ------------------------------------------------------------------------------
# Evaluation Software for NEMD simulations of same state - EvalStateNEMD
# Functions to calculate Newtonian viscosity for state
# ---
# created by Jens Wagner, 29.10.2021
# ------------------------------------------------------------------------------

# Main Function
function eval_state_shear(subfolder::Array{String,1}, inpar::Opts)

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
    T, p, ρ, x = StaticProperties(subfolder)

    # Get NEMD results of single simulations
    η_all, s_rate_all, T_all, p_all, ρ_all = collectNEMD(subfolder)
    η = getfield.(η_all,:val)
    s_rate = getfield.(s_rate_all,:val)

    # Define functions#
    fun_Carreau(x,p) = p[1] .* (1 .+ (x ./ p[2]).^2).^((p[3].-1)./2)    # p[1]: η_N, p[2]: s_rate_0, p[3]: n
    fun_Carreau_log(x,p) = log.(fun_Carreau(x,p))
    fun_Eyring(x,p) = p[2] ./ x .* asinh.(p[1] .* x ./ p[2])            # p[1]: η_N, p[2]: σ_E
    fun_Eyring_log(x,p) = log.(fun_Eyring(x,p))

    s_rate_v = 10 .^ [LinRange(log10(minimum(s_rate)),log10(maximum(s_rate)),100);]

    # Fit Carreau
    fit_Carreau = curve_fit(fun_Carreau_log, s_rate, log.(η), [η[findmin(s_rate)[2]], median(s_rate),0.5], x_tol=1e-20, g_tol=1e-20)
    SSE_Carreau = sum( (log.(η) .- fun_Carreau_log(s_rate,fit_Carreau.param)).^2 )

    # Fit Eyring
    fit_Eyring = curve_fit(fun_Eyring, s_rate, η, [η[findmin(s_rate)[2]], median(s_rate.*η)], x_tol=1e-20, g_tol=1e-20)
    SSE_Eyring = sum( (log.(η) .- fun_Eyring_log(s_rate,fit_Eyring.param)).^2 )

    # Figure γ-η
    figure()
    errorbar(s_rate, η, fmt="ok", mfc=:white, yerr=getfield.(η_all,:err), label="simulations", ecolor="k", elinewidth = 1, capsize = 2)
    plot(s_rate_v, fun_Carreau(s_rate_v,fit_Carreau.param), "b--", label="Carreau fit")
    plot(s_rate_v, fun_Eyring(s_rate_v,fit_Eyring.param), "r:", label="Eyring fit")
    if !(reduced_units)
        title("Viscosity \$\\eta(\\dot{\\gamma})\$ (\$T = $(round(T.val, digits=2)){\\rm~K}\$, \$\\rho = $(round(ρ.val,digits=5)){\\rm~g~ml^{-1}}\$)")
        xlabel("\$\\dot{\\gamma}{\\rm~/~s^{-1}}\$")
        ylabel("\$\\eta{\\rm~/~Pa~s}\$")
    elseif reduced_units
        title("Viscosity \$η^{*}(\\dot{\\gamma}^{*})\$ (\$T^{*} = $(round(T.val, digits=2))\$, \$\\rho^{*} = $(round(ρ.val,digits=5))\$")
        xlabel("\$\\dot{\\gamma}^{*}\$")
        ylabel("\$\\eta^{*}\$")
    end
    xscale("log");  yscale("log")
    tick_params(which="both", right=true, top=true ,direction="in")
    legend();       tight_layout()
    savefig(string(folder,"/fig_shearrate_viscosity.pdf"))
    close()

    # Choose the fit with the lower SSE
    if SSE_Carreau < SSE_Eyring
        best_fit = "Carreau"
        η_N = SingleDat(fit_Carreau.param[1], stderror(fit_Carreau)[1], margin_error(fit_Carreau,0.05)[1])
    else
        best_fit = "Eyring"
        η_N = SingleDat(fit_Eyring.param[1], stderror(fit_Eyring)[1], margin_error(fit_Eyring,0.05)[1])
    end
    note_bestfit = " (Newtonian viscosity calculated by the $best_fit model)"

    # Create result
    res = ResultsDatNEMD(T, p, ρ, x, [], [], [], η_N, [], [])
    output_resultsNEMD(res, folder; note=note_bestfit)

    create_parameterfile_shear(folder, fit_Carreau.param, fit_Eyring.param)

    fID = open("$folder/sim_dat.csv","w")
    println(fID,"T,p,rho,eta,s_rate")
    writedlm(fID,[getfield.(T_all,:val) getfield.(p_all,:val) getfield.(ρ_all,:val) η s_rate],",")
    close(fID)
end

## Subfunctions ----------------------------------------------------------------
# Collect NEMD properties from single simulations
function collectNEMD(subf)
    n = length(subf)
    T_all = Array{SingleDat,1}(undef,n)
    p_all = Array{SingleDat,1}(undef,n)
    ρ_all = Array{SingleDat,1}(undef,n)
    η_all = Array{SingleDat,1}(undef,n)
    s_rate_all = Array{SingleDat,1}(undef,n)

    for i = 1:n
        res = load_result_NEMD(string(subf[i],"/result.dat"))
        T_all[i] = res.T
        p_all[i] = res.p
        ρ_all[i] = res.ρ
        η_all[i] = res.η
        s_rate_all[i] = res.s_rate
    end

    return η_all, s_rate_all, T_all, p_all, ρ_all
end

# Create fit parameter file
function create_parameterfile_shear(folder, p_Carreau, p_Eyring)

    # Create Filepath
    path = string(folder,"/fit_parameters.dat")

    # Write to file
    fID = open(path,"w")
    line1 = "# Created by MD - Bulk Evaluation, folder: $folder, time: $(Dates.format(now(),"yyyy-mm-dd HH:MM:SS"))"
    line2 = "# Parameters of fit for Carreau and Eyring model"
    if !(reduced_units)
        line3 = "# units: [η_N] = Pa*s, [γ_0] = 1/s, [n] = 1, [σ_E] = Pa"
    else
        line3 = "# units: reduced"
    end
    header = string(line1,"\n",line2,"\n",line3,"\n")
    print(fID,header)

    # Write data
    line3 = "# ---------- Carreau model ---------- \n"
    print(fID,line3)
    print_prop(fID, p_Carreau[1], "η_N")
    print_prop(fID, p_Carreau[2], "γ_0")
    print_prop(fID, p_Carreau[3], "n")
    line4 = "# ---------- Eyring model ---------- \n"
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
    res = ResultsDatNEMD([],[],[],[],[],[],[],[],[],[])

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
                if (name == "T")        res.T = SingleDat(val,std,err) end
                if (name == "p")        res.p = SingleDat(val,std,err) end
                if (name == "ρ")        res.ρ = SingleDat(val,std,err) end
                if (name == string("x",num))
                                        res.x = vcat(res.x, SingleDat(val,std,err))
                                        if length(res.x) != num @warn("res.x doesn't fit num!") end
                end
                if (name == "Etot")     res.Etot = SingleDat(val,std,err) end
                if (name == "Ekin")     res.Ekin = SingleDat(val,std,err) end
                if (name == "Epot")     res.Epot = SingleDat(val,std,err) end
                if (name == "η")        res.η = SingleDat(val,std,err) end
                if (name == "s_rate")   res.s_rate = SingleDat(val,std,err) end
                if (name == "λ")        res.λ = SingleDat(val,std,err) end
            end
        end
    end
    return res
end

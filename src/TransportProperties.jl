## TransportProperties.jl
# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - EvalData
# Function to apply Time Decomposition Method to Calculate Viscosity
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------
# References:
# [1] Fischer, M.; Bauer, G.; Gross, J. Force Fields with Fixed Bond Lengths and with Flexible Bond Lengths: Comparing Static and Dynamic Fluid Properties. Journal of Chemical & Engineering Data 2020. https://doi.org/10.1021/acs.jced.9b01031.
# [2] Maginn, E. J.; Messerly, R. A.; Carlson, D. J.; Roe, D. R.; Elliott, J. R. Best Practices for Computing Transport Properties 1. Self-Diffusivity and Viscosity from Equilibrium Molecular Dynamics [Article v1.0]. LiveCoMS 2019, 1 (1). https://doi.org/10.33011/livecoms.1.1.6324.

## Main
function TransportProperties(T::Float64,L_box::Float64,set::set_TDM)
    # Set what to calculate and output mode
    do_η = 1
    do_D = 1
    do_λ = 1
    out_mode = 1        # 1 - Take mean value of bootstrapping
                        # 2 - Take result of single fit

    # Create new folder if it not yet exists
    if !(isdir(string(set.folder,"/TransportProperties")))
        mkdir(string(set.folder,"/TransportProperties"))
    end

    fID_info = open(string(set.folder,"/TransportProperties/info.txt"),"w")
    println(fID_info,"Start: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS"))

    # Load existing result file
    do_calc = 1                 # 1 - Calculation even if it already exists
                                # 0 - Take old result if exists
    if isfile(string(set.folder,"/result.dat"))
        res = load_result(string(set.folder,"/result.dat"))
    else
        do_calc = 1
    end

    # VISCOSITY
    if do_η == 1 && (do_calc == 1 || !(typeof(res.η) == single_dat))
        # Load runs
        ηmat, t = load_runs("viscosity.dat", 2, set)
        # Calculation of viscosity value by TDM method
        set_η = set
        set_η.name = "Viscosity"
        set_η.symbol = "\\eta"
        if !(reduced_units) set_η.unit = "(Pa*s)" elseif (reduced_units) set_η.unit = "-" end
        set_η.do_out = true
        plot_all_curves(t, ηmat, set_η)
        ηmat = eliminate_outlier_curves(ηmat,set_η,fID_info)
        ηfit, set_η = TDM(ηmat, t, set_η)
        # Calculation of statistical uncertainties by bootstrapping method
        if set.nboot > 0
            println("Bootstrapping Viscosity ...")
            ηboot, ηstd, ηerr = bootstrapping(ηmat, t, set)
        else
            ηboot = NaN
            ηstd = NaN
            ηerr = NaN
        end
        # Save in struct
        if !(isnan(ηboot)) && out_mode == 1
            η = single_dat(ηboot, ηstd, ηerr)
        else
            η = single_dat(ηfit, ηstd, ηerr)
        end

    elseif do_η == 0
        η = []
    else
        η = res.η
    end

    println(fID_info,"Viscosity DONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS"))

    # DIFFUSION COEFFICIENT
    if do_D == 1 && (do_calc == 1 || !(typeof(res.D) == single_dat))
        D = Array{single_dat,1}(undef,0)
        imol = 0
        Nmol = 10
        while true
            imol += 1
            Dv = Float64[]
            for i = 1:length(set.subfolder)
                file = string(set.subfolder[i],"/result.dat")
                if isfile(file)
                    restmp = load_result(file)
                    append!(Dv,restmp.D[imol].val)
                    Nmol = length(restmp.x)
                end
            end
            # Finite size correction [1]
            ξ = 2.837298
            if (reduced_units)
                Dcorr = T*ξ/(6*π*η.val*L_box)
            else
                Dcorr = kB*T*ξ/(6*π*η.val*L_box)
            end
            Dval = mean(Dv) + Dcorr
            D = vcat(D,single_dat(Dval,std(Dv),std(Dv)/sqrt(length(Dv))))
            if imol == Nmol break end
        end
    elseif do_D == 0
        D = []
    else
        D = res.D
    end

    println(fID_info,"Self-Diffusion Coef. DONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS"))

    # THERMAL CONDUCTIVITY
    if do_λ == 1 && (do_calc == 1 || !(typeof(res.λ) == single_dat))
        # Load runs
        λmat, t = load_runs("thermalconductivity.dat", 2, set)
        # Calculation of viscosity value by TDM method
        set_λ = set
        set_λ.name = "Thermal Conductivity"
        set_λ.symbol = "\\lambda"
        if !(reduced_units) set_λ.unit = "/(W/(m*K))" elseif (reduced_units) set_λ.unit = "-" end
        set_λ.do_out = true
        plot_all_curves(t, λmat, set_λ)
        λmat = eliminate_outlier_curves(λmat,set_λ,fID_info)
        λfit, set_λ = TDM(λmat, t, set_λ)
        # Calculation of statistical uncertainties by bootstrapping method
        if set.nboot > 0
            println("Bootstrapping Thermal Conductivity ...")
            λboot, λstd, λerr = bootstrapping(λmat, t, set)
        else
            λboot = NaN
            λstd = NaN
            λerr = NaN
        end
        # Save in struct
        if !(isnan(λboot)) && out_mode == 1
            λ = single_dat(λboot, λstd, λerr)
        else
            λ = single_dat(λfit, λstd, λerr)
        end
    elseif do_λ == 0
        λ = []
    else
        λ = res.λ
    end

    println(fID_info,"Thermal Conductivity DONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS"))
    println(fID_info,"End: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS"))
    close(fID_info)

    η_V = []
    return η, η_V, D, λ
end

# Subfunctions
# Load run data
function load_runs(file::String, index::Int64, set::set_TDM)
    n = length(set.subfolder)
    tdat = Array{Array{Float64,1}}(undef,n)
    dat = Array{Array{Float64,1}}(undef,n)
    for i = 1:n
        filedat = readdlm(string(set.subfolder[i],"/",file), skipstart=3)
        dat[i] = filedat[:,index]
        tdat[i] = filedat[:,1]
    end
    L = minimum(length.(tdat))
    tdat = getindex.(tdat, repeat([1:L],n))
    t = unique(tdat)[1]
    if (size(unique(tdat),2) > 1) error("Data don't match!") end
    mat = hcat(getindex.(dat, repeat([1:L],n))...)

    return mat, t
end

# Function to perform 1 TDM calculation
@everywhere function TDM(mat::Array{Float64,2}, t::Array{Float64,1}, set::set_TDM)
    t = t .- t[1]
    n = size(mat,2)
    ave_t = mean(mat, dims=2)
    std_t = sqrt.( sum((mat .- ave_t).^2, dims=2) ./(n-1))[:,1]

    # Function to fit running standard deviation (p[1]: A, p[2]: b)
    @. fun_std(x, p) = p[1] .* x .^ p[2]
    # Function to fit running viscosity (p[1]: η(t→∞), p[2]: α, p[3]: β₁, p[4]: β₂), following Ref. (1)
    @. fun_ave(x, p) = p[1] .* ( ((p[2].*p[3].*(1 .-exp(-x./p[3])) .+
                                  (1 .-p[1]).*p[4]).*(1 .-exp(-x./p[4]))) ./
                                  (p[2].*p[3].+(1 .-p[1]).*p[4]) )

    # Fit standard deviation
    k = 0
    converged = false
    p0_std = [[1.0,1.0], [2.0,0.5], [0.5,2.0], [1e-4,1.0], [1.0,1e-4]]
    fit_std = []
    pos = findfirst((std_t./ave_t)[:] .> set.cutcrit)
    if isnothing(pos) pos = length(t) end
    cut_std = minimum([length(t), round(Int64,1.1 .* pos)])
    while !(converged)
        k += 1
        fit_std = curve_fit(fun_std, t[1:cut_std], std_t[1:cut_std], p0_std[k])
        converged = fit_std.converged
        if (k == length(p0_std)) break end
    end
    if !(fit_std.converged) println("Std: Not converged!") end

    # Calculation of tcut (or cut)
    skip_std = findfirst(t .>= 2)
    cut = findfirst(fun_std(t[skip_std:end],fit_std.param)./ave_t[skip_std:end] .> set.cutcrit)
    if isnan(set.tskip)
        skip = round(Int64,cut*0.02)
    else
        skip = findfirst(t .>= set.tskip)
    end
    if fit_std.converged && !(isnothing(cut))
        cut += skip-1
        set.tcut = t[cut]

        # Fit η_ave(t) (start from tstart ps)
        k_exponent = 1.0        # 1 for original TDM; if bad fit -> reduce to ~ 0.1 - 1.0
        w = 1 ./ t[skip:cut].^(fit_std.param[2]*k_exponent)
        k = 0
        k_ex = 0
        converged = false

        p0_1 = mean(ave_t[round(Int64,cut/2):cut])
        p0_ave = [  [p0_1, 1.0, 1.0,  1.0],
                    [p0_1, 1.0, 0.1,  10.0],
                    [p0_1, 1.0, 1e-3, 1.0],
                    [p0_1, 0.1, 0.1,  1.0],
                    [0.1,  1.0, 0.1,  1.0] ]
        fit_ave = []

        while !(converged)
            k += 1
            try
                fit_ave = curve_fit(fun_ave, t[skip:cut], ave_t[skip:cut], w, p0_ave[k])
                converged = fit_ave.converged
            catch x
                if  string(x) == "LinearAlgebra.SingularException(4)" && k_ex < 5
                    for i=1:length(p0_ave) setindex!(p0_ave,p0_ave[i].*[2,1,1,1],i) end
                    k -= 1;         k_ex += 1
                else rethrow(x)
                end
            end
            if (k == length(p0_ave)) break end
        end

        val = fit_ave.param[1]

        if !(fit_ave.converged)
            val = NaN
            println("Ave: Not converged!")
        end
    else
        val = NaN
    end

    # Output plot
    if set.do_out && !isnan(val)
        outfolder = string(set.folder,"/TransportProperties/")

        # Save all single simulation data and the averaged data of the transport property
        line1 = string("# Created by MD - Bulk Evaluation, Folder: ", set.folder)
        line2 = string("# t[ps] ave[",set.unit,"]")
        if (reduced_units) line2 = string("# t* ave*") end
        for i = 1:size(mat,2) line2 = string(line2," sim",i,"[",set.unit,"]") end
        header = string(line1,"\n",line2)
        file = string(outfolder,set.name,".dat")
        fID = open(file,"w"); println(fID,header)
        writedlm(fID, hcat(t[1:cut],ave_t[1:cut],mat[1:cut,:]), " ")
        close(fID)

        # Plot: γ(t) average values and fitted curve
        figure()
        if !(reduced_units)
            xlabel(L"t / ps")
            ylabel(latexstring(set.symbol," / ",set.unit))
        elseif reduced_units
            xlabel(L"t*")
            ylabel(latexstring(set.symbol,"*"))
        end
        plot(t[1:cut],mat[1:cut,:], color="gray", linestyle=":", linewidth=0.1)
        plot(t[1:cut],ave_t[1:cut], color="blue", label="Average", linewidth=0.1)
        plot(t[1:cut],fun_ave(t[1:cut],fit_ave.param), color="red", label="Fit", linewidth=1)
        title(string(latexstring("Fit: ",set.symbol,"_0 = ",fit_ave.param[1],", \\alpha = ",fit_ave.param[2]),",\n",
                     latexstring("\\beta_1 = ",fit_ave.param[3],", \\beta_2 = ",fit_ave.param[4])),fontsize=10)
        legend(loc="upper left")
        tight_layout()
        savefig(string(outfolder,set.name,".pdf"))

        # Plot: σ_γ(t) standard deviation and fitted curve
        figure()
        if !(reduced_units)
            xlabel(L"t / ps")
            ylabel(latexstring("\\sigma_{",set.symbol,"} / ",set.unit))
        elseif reduced_units
            xlabel(L"t*")
            ylabel(latexstring("\\sigma_{",set.symbol,"}*"))
        end
        plot(t[1:cut],std_t[1:cut], color="blue", label="Average", linewidth=0.1)
        plot(t[1:cut],fun_std(t[1:cut],fit_std.param), color="red", label="Fit", linewidth=1)
        title(string(latexstring("Fit: A = ",fit_std.param[1],", b = ",fit_std.param[2])),fontsize=10)
        legend(loc="upper left")
        tight_layout()
        savefig(string(outfolder,set.name,"_std.pdf"))

        if (val < 0)
            error("Negative value for ",set.name,"!")
        elseif isnan(val)
            error("NaN value for ",set.name," (-> fit did not converge)!")
        end
    end
    close("all")

    return val, set
end

# Function to apply bootstrpaping on data
function bootstrapping(mat, t, set)
    nboot = set.nboot
    nsim = size(mat,2)

    # Create matrix of #nboot random combinations
    randmat = Array{Array{Int64,1},1}(undef,Int64(2*nboot))
    for i = 1:Int64(2*nboot)
        k = 0
        while k < 10
            k += 1
            randmat[i] = sort(rand(1:nsim,nsim))
            if (length(unique(randmat[i])) > nsim/3) break end
        end
    end
    bootmat = unique(randmat)[1:nboot,:]

    # Create cut off for every TDM evaluation
    cutcrit = rand(Float64, nboot) .* 0.2 .+ set.cutcrit

    # TDM for all combinations of single simulations
    vals = pmap((x1,x2)->TDMboot(mat,t,set,x1,x2), bootmat, cutcrit)

    # Remove outliers
    k = 3
    what_NaN = abs.(vals) .> k*median(vals[(!).(isnan.(vals))])
    vals[what_NaN] = NaN .* ones(sum(what_NaN))
    what_NaN = abs.(vals) .< median(vals[(!).(isnan.(vals))])/k
    vals[what_NaN] = NaN .* ones(sum(what_NaN))

    # # Calculate error from distribution of TDM values
    # nb = round(Int,nboot/20)            # number of bins
    # h = fit(Histogram, vals[.!isnan.(vals)], nbins = nb)
    # edg = Array(h.edges[1])
    # dedg = mean(edg[2:end] .- edg[1:end-1])
    # pts = (edg[2:end].+edg[1:end-1])./2
    # wgs = h.weights ./ (sum(h.weights) * dedg)

    fit_normal = fit_mle(Normal, vals[.!isnan.(vals)])
    μ = fit_normal.μ
    σ = fit_normal.σ
    x = Array(LinRange(minimum(vals),maximum(vals),100))
    y = exp.(-(x.-μ).^2 ./ (2 .*σ.^2)) ./ sqrt.(2 .*π.*σ.^2)

    # Output to dat file
    fID = open(string(set.folder,"/TransportProperties/bootstrapping_",set.name,".dat"),"w")
    writedlm(fID,vals," ")
    close(fID)

    # Plot histogram
    figure()
    hist(vals, bins=round(Int,nboot/15), label="Data", density=true)
    plot(x,y, label="Fit", linewidth=2)
    outfolder = string(set.folder,"/TransportProperties/")
    legend(loc="upper left")
    title(string(L"Histogram with: $\mu = $",μ," ",set.unit,L", $\sigma = $",σ," ",set.unit))
    tight_layout()
    savefig(string(outfolder,"histogram_",set.name,".pdf"))
    close()

    err = 1.96*σ

    return μ, σ, err
end

# Help function for efficient TDM calculation in bootstrapping
@everywhere function TDMboot(mat,t,set,bootsample,crit)
    set.cutcrit = crit
    set.do_out = false
    matTMP = mat[:,bootsample]
    val, set = TDM(matTMP,t,set)
    return val
end

# Function to eliminate outlier curves
function eliminate_outlier_curves(mat::Array{Float64,2},set::set_TDM,fID)
    # Settings
    fraction = 0.25      # Fraction what timesteps to use for mean calculation
    factor = 3          # Factor to define outlier curves

    nrow = size(mat,1)
    ncol = size(mat,2)

    # Calculation of means
    means = mean(mat[1:Int(round(fraction*nrow)),:],dims=1)
    # Calculation of median of means
    m = median(means[(!).(isnan.(means))])

    # Reduction of matrix
    index = (abs.(means) .< (factor*m))[:]
    mat = mat[:,index]

    if sum(index .== 0) > 0
        println(set.name,": Elimination of ",sum(index .== 0)," simulation(s) - Sims.: ",findall(index .== 0))
        println(fID,set.name,": Elimination of ",sum(index .== 0)," simulation(s) - Sims.: ",findall(index .== 0))
    end
    return mat
end

# Function to plot all curves to a figure
function plot_all_curves(t, mat, set)
    # Output folder
    outfile = string(set.folder,"/TransportProperties/",set.name,"_all.pdf")

    figure()
    for i = 1:size(mat,2)
        plot(t,mat[:,i],label=string("Sim ",i),lw=0.5)
    end
    if !(reduced_units)
        xlabel(L"t / ps")
        ylabel(latexstring(set.symbol," / ",set.unit))
    elseif reduced_units
        xlabel(L"t*")
        ylabel(latexstring(set.symbol,"*"))
    end
    legend(loc="upper left", fontsize=5, ncol=maximum([1,round(Int,size(mat,2)/10)]))
    tight_layout()
    savefig(outfile)
    close()
end

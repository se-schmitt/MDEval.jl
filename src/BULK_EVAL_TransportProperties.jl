# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - EvalData
# Function to apply Time Decomposition Method to Calculate Viscosity
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

# Main
function TransportProperties(state::state_info,set::set_TDM)
    # Set what to calculate
    do_η = 1
    do_D = 1
    do_λ = 1

    if set.nboot == 0
        do_η = 0;   do_D = 0;   do_λ = 0
    end

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
        ηval, set_η = TDM(ηmat, t, set_η)
        # Calculation of statistical uncertainties by bootstrapping method
        if set.nboot > 0
            println("Bootstrapping Viscosity ...")
            ηstd, ηerr = bootstrapping(ηmat, t, set)
        else
            ηstd = NaN
            ηerr = NaN
        end
        # Save in struct
        η = single_dat(ηval, ηstd, ηerr)
    elseif do_η == 0
        η = []
    else
        η = res.η
    end

    println(fID_info,"Viscosity DONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS"))

    # DIFFUSION COEFFICIENT
    if do_D == 1 && (do_calc == 1 || !(typeof(res.D) == single_dat))
        Dv = []
        for i = 1:length(set.subfolder)
            file = string(set.subfolder[i],"/result.dat")
            if isfile(file)
                restmp = load_result(file)
                append!(Dv,restmp.D.val)
            end
        end
        ξ = 2.837298
        L = (state.m /state.ρ * 1e-6)^(1/3)
        Dcorr = kB*state.T*ξ/(6*π*η.val*L)
        Dval = mean(Dv) .+ Dcorr
        D = single_dat(Dval,std(Dv),std(Dv)/sqrt(length(Dv)))
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
        λval, set_λ = TDM(λmat, t, set_λ)
        # Calculation of statistical uncertainties by bootstrapping method
        if set.nboot > 0
            println("Bootstrapping Thermal Conductivity ...")
            λstd, λerr = bootstrapping(λmat, t, set)
        else
            λstd = NaN
            λerr = NaN
        end
        # Save in struct
        λ = single_dat(λval, λstd, λerr)
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
        filedat = readdlm(string(set.subfolder[i],"/",file), skipstart=2)
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
    # Function to fit running viscosity (p[1]: η(t→∞), p[2]: α, p[3]: β₁, p[4]: β₂)
    @. fun_ave(x, p) = p[1] .* ( ((p[2].*p[3].*(1 .-exp(-x./p[3])) .+
                                  (1 .-p[1]).*p[4]).*(1 .-exp(-x./p[4]))) ./
                                  (p[2].*p[3].+(1 .-p[1]).*p[4]) )

    skip = findfirst(t .>= set.tskip)

    # Fit standard deviation
    k = 0
    converged = false
    p0_std = [[1.0,1.0], [2.0,0.5], [0.5,2.0], [1e-4,1.0], [1.0,1e-4]]
    fit_std = []
    while !(converged)
        k += 1
        fit_std = curve_fit(fun_std, t, std_t, p0_std[k])
        converged = fit_std.converged
        if (k == length(p0_std)) break end
    end
    if !(fit_std.converged) println("Std: Not converged!") end

    # Calculation of tcut (or cut)
    cut = findfirst(fun_std(t[skip:end],fit_std.param)./ave_t[skip:end] .> set.cutcrit)
    if fit_std.converged && !(isnothing(cut))
        cut += skip-1
        set.tcut = t[cut]

        # Fit η_ave(t) (start from tstart ps)
        w = 1 ./ t[skip:cut].^fit_std.param[2]
        k = 0
        converged = false
        p0_1 = mean(ave_t[round(Int64,cut/2):cut])
        p0_ave = [  [p0_1, 1.0, 1.0,  1.0],
                    [p0_1, 1.0, 0.1,  10.0],
                    [p0_1, 1.0, 1e-3, 1e3],
                    [p0_1, 0.1, 0.1, 1.0],
                    [0.1,  1.0, 0.1, 1.0] ]
        fit_ave = []
        while !(converged)
            k += 1
            fit_ave = curve_fit(fun_ave, t[skip:cut], ave_t[skip:cut], w, p0_ave[k])
            converged = fit_ave.converged
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
        tight_layout()
        if !(reduced_units) xlabel(L"t / ps") elseif reduced_units xlabel(L"t*") end
        nstring = latexstring(set.symbol," / ",set.unit)
        rstring = latexstring(set.symbol,"*")
        if !(reduced_units) ylabel(nstring) elseif reduced_units ylabel(rstring) end
        plot(t[1:cut],mat[1:cut,:], color="gray", linestyle=":", linewidth=0.1)
        plot(t[1:cut],ave_t[1:cut], color="blue", label="Average", linewidth=0.1)
        plot(t[1:cut],fun_ave(t[1:cut],fit_ave.param), color="red", label="Fit", linewidth=1) #label="Fit",linecolor=:red
        legend(loc="upper left")
        splot = string(outfolder,set.name,".pdf")
        savefig(splot)
        close()

        if (val < 0)
            error("Negative value for ",set.name,"!")
        elseif isnan(val)
            error("NaN value for ",set.name," (-> fit did not converge)!")
        end
    end

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
    cutcrit = rand(Float64, nboot) .* 0.2 .+ 0.4

    # TDM for all combinations of single simulations
    vals = pmap((x1,x2)->TDMboot(mat,t,set,x1,x2), bootmat, cutcrit)

    # Remove outliers
    what_NaN = abs.(vals) .> 2*median(vals)
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
    tight_layout()
    hist(vals, bins=round(Int,nboot/15), label="Data")
    plot(x,y, label="Fit", linewidth=2)
    outfolder = string(set.folder,"/TransportProperties/")
    splot = string(outfolder,"histogram_",set.name,".pdf")
    legend(loc="upper left")
    savefig(splot)
    close()

    err = 1.96*σ

    return σ, err
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
    fraction = 0.5      # Fraction what timesteps to use for mean calculation
    factor = 2          # Factor to define outlier curves

    nrow = size(mat,1)
    ncol = size(mat,2)

    # Calculation of means
    means = mean(mat[1:Int(round(fraction*nrow)),:],dims=1)
    # Calculation of median of means
    m = median(means)

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
    if !(reduced_units) xlabel(L"t / ps") elseif reduced_units xlabel(L"t*") end
    nstring = latexstring(set.symbol," / ",set.unit)
    rstring = latexstring(set.symbol,"*")
    if !(reduced_units) ylabel(nstring) elseif reduced_units ylabel(rstring) end
    legend(loc="upper left", fontsize=5, ncol=round(Int,size(mat,2)/10))
    tight_layout()
    savefig(outfile)
    close()
end

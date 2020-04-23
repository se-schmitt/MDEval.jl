# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - EvalData
# Function to apply Time Decomposition Method to Calculate Viscosity
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

# Main
function TransportProperties(set::set_TDM)
    if !(isdir(string(set.folder,"/TransportProperties")))
        mkdir(string(set.folder,"/TransportProperties"))
    end

    # VISCOSITY
    # Load runs
    ηmat, t = load_runs("viscosity.dat", 2, set)
    # Calculation of viscosity value by TDM method
    set_η = set
    set_η.name = "Viscosity"
    set_η.unit = "Pa*s"
    set_η.do_out = true
    ηval, set_η = TDM(ηmat, t, set_η)
    # Calculation of statistical uncertainties by bootstrapping method
    println("Bootstrapping Viscosity ...")
    ηstd, ηerr = bootstrapping(ηmat, t, set)
    # Save in struct
    η = single_dat(ηval, ηstd, ηerr)

    # DIFFUSION COEFFICIENT
    D = single_dat([],[],[])

    # THERMAL CONDUCTIVITY
    # Load runs
    λmat, t = load_runs("thermalconductivity.dat", 2, set)
    # Calculation of viscosity value by TDM method
    set_λ = set
    set_λ.name = "Thermal Conductivity"
    set_λ.unit = "W/(m*K)"
    set_λ.do_out = true
    λval, set_λ = TDM(λmat, t, set_λ)
    # Calculation of statistical uncertainties by bootstrapping method
    println("Bootstrapping Thermal Conductivity ...")
    λstd, λerr = bootstrapping(λmat, t, set)
    # Save in struct
    λ = single_dat(λval, λstd, λerr)

    η_V = single_dat([],[],[])
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
    fit_std = curve_fit(fun_std, t, std_t, [1.0,1.0])
    if !(fit_std.converged) println("Std: Not converged!") end
    # Calculation of tcut (or cut)
    cut = findfirst(fun_std(t[skip:end],fit_std.param)./ave_t[skip:end] .> set.cutcrit)
    if fit_std.converged && !(isnothing(cut))
        cut += skip-1
        set.tcut = t[cut]

        # Fit η_ave(t) (start from tstart ps)
        w = 1 ./ t[skip:cut].^fit_std.param[2]
        p0 = [ave_t[cut],1.0,1.0,10.0]
        fit_ave = curve_fit(fun_ave, t[skip:cut], ave_t[skip:cut], w, p0)

        val = fit_ave.param[1]

        # println("cut_crit: ",set.cutcrit,"; cut pos: ",cut,"; t_cut: ",set.tcut,
        #         "; fit ave: ",fit_ave.converged,"; val: ",val)
        if !(fit_ave.converged)
            val = NaN
        end
    else
        val = NaN
    end

    # Output plot
    if set.do_out
        outfolder = string(set.folder,"/TransportProperties/")
        valstr = string(round(val/10^floor(log10(val)),digits=3),"e",Int64(floor(log10(val))))
        plt = plot(t[1:cut],[ave_t[1:cut],fun_ave(t[1:cut],fit_ave.param)],
            xlabel="t / ps",ylabel=string(set.name," / ",set.unit), dpi=400, legend=:bottomright,
            label=[string("Averaged ",set.name) string("Fit: ",set.name," = ",valstr," ",set.unit)])
        png(plt,string(outfolder,set.name,".png"))
        # fID = open(string(outfolder,set.name,".info"),"w")
        # println("Folder: ",set.folder)
    end

    return val, set
end

# Function to apply bootstrpaping on data
function bootstrapping(mat, t, set)
    nboot = set.nboot
    nsim = size(mat,2)

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

    cutcrit = rand(Float64, nboot) .* 0.2 .+ 0.4

    vals = pmap((x1,x2)->TDMboot(mat,t,set,x1,x2), bootmat, cutcrit)

    nb = 50 # number of bins
    h = fit(Histogram, vals[.!isnan.(vals)], nbins = nb)
    edg = Array(h.edges[1])
    dedg = mean(edg[2:end] .- edg[1:end-1])
    pts = (edg[2:end].+edg[1:end-1])./2
    wgs = h.weights ./ (sum(h.weights) * dedg)

    fit_normal = fit_mle(Normal, vals[.!isnan.(vals)])
    μ = fit_normal.μ
    σ = fit_normal.σ
    x = Array(LinRange(minimum(pts),maximum(pts),100))
    y = exp.(-(x.-μ).^2 ./ (2 .*σ.^2)) ./ sqrt.(2 .*π.*σ.^2)

    plt = histogram(vals,normalize=true,dpi=400,label="Data")
    plot!(x,y,label="Fit")
    outfolder = string(set.folder,"/TransportProperties/")
    png(string(outfolder,"histogram_",set.name,".png"))

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

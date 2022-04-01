## GreenKubo.jl
# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - GreenKubo
# Containing functions to apply Green-Kubo method to viscosity and thermal
# conducitvity
# ---
# created by Sebastian Schmitt, 13.03.2021
# ------------------------------------------------------------------------------
# References:
# (1) Humbert, M. T.; Zhang, Y.; Maginn, E. J. PyLAT: Python LAMMPS Analysis Tools. J. Chem. Inf. Model. 2019, 59 (4), 1301–1305. https://doi.org/10.1021/acs.jcim.9b00066.
# (2) Maginn, E. J.; Messerly, R. A.; Carlson, D. J.; Roe, D. R.; Elliott, J. R. Best Practices for Computing Transport Properties 1. Self-Diffusivity and Viscosity from Equilibrium Molecular Dynamics [Article v1.0]. LiveCoMS 2019, 1 (1). https://doi.org/10.33011/livecoms.1.1.6324.
# (3) Fischer, M.; Bauer, G.; Gross, J. Force Fields with Fixed Bond Lengths and with Flexible Bond Lengths: Comparing Static and Dynamic Fluid Properties. Journal of Chemical & Engineering Data 2020. https://doi.org/10.1021/acs.jced.9b01031.

## Viscosity -------------------------------------------------------------------
# Calculation of viscosity from pressure tensor
function calc_viscosities(info::info_struct, mode::String; mode_acf::String, CorrLength::Int64=0, SpanCorrFun::Int64=1, nEvery::Int64=1)
    # Loading Pressure File
    dat = load_pressure(info, nEvery)
    if (nEvery > 1) print("\n"); @warn("Only every $(nEvery). step of pressure tensor used for ACF calculation.") end

    if (dat isa pressure_dat)
        what = (dat.step .>= info.n_equ)

        # Unit conversion and GK factor
        metal2Pas = 1e-32
        factor = metal2Pas .* mean(dat.V[what]) ./ (kB .* mean(dat.T[what]))

        # Calculation of T and ρ
        T_std_err = block_average(dat.T[what],N_blocks=info.n_blocks)
        T = single_dat(mean(dat.T[what]), T_std_err[1], T_std_err[2])
        ρ_std_err = block_average(info.natoms ./ dat.V[what],N_blocks=info.n_blocks)
        ρ = single_dat(mean(info.natoms ./ dat.V[what]), ρ_std_err[1], ρ_std_err[2])

        if (reduced_units) factor = mean(dat.V[what]) ./ mean(dat.T[what]) end

        # Define start and end steps of individual correlation functions
        pos_step_start = []
        pos_step_end = []
        if mode == "single"
            for i_step in info.n_equ:SpanCorrFun:(maximum(dat.step) - CorrLength)
                append!(pos_step_start, findfirst(dat.step .>= i_step))
                append!(pos_step_end,   findfirst(dat.step .>= i_step+CorrLength))
            end
        elseif mode == "tdm"
            append!(pos_step_start,findfirst(what))
            append!(pos_step_end,length(dat.step))
        end

        # Allocate arrays
        η_t_all =     NaN .* ones(pos_step_end[1] .- pos_step_start[1]+1, length(pos_step_start))
        η_V_t_all =   NaN .* ones(pos_step_end[1] .- pos_step_start[1]+1, length(pos_step_start))
        acf_η_all =   NaN .* ones(pos_step_end[1] .- pos_step_start[1]+1, length(pos_step_start))
        acf_η_V_all = NaN .* ones(pos_step_end[1] .- pos_step_start[1]+1, length(pos_step_start))

        for i in 1:length(pos_step_start)
            what_pos = pos_step_start[i]:pos_step_end[i]

            # Shear viscosity
            # Calculate autocorrelation function
            P_mat = [dat.pxy[what_pos],
                     dat.pxz[what_pos],
                     dat.pyz[what_pos],
                     dat.pxx[what_pos] .- dat.pyy[what_pos],
                     dat.pyy[what_pos] .- dat.pzz[what_pos],
                     dat.pxx[what_pos] .- dat.pzz[what_pos]]
            P_diag = [dat.pxx[what_pos], dat.pyy[what_pos], dat.pzz[what_pos]]

            acf_P = pmap(x -> calc_acf(x, mode_acf), vcat(P_mat,P_diag))

            # Shear viscosity
            # -> method adopted from PyLat code (1)
            acf_η_tmp = (acf_P[1] .+ acf_P[2] .+ acf_P[3])./6 .+
                        (acf_P[4] .+ acf_P[5] .+ acf_P[6])./24
            acf_η_all[:,i] = acf_η_tmp
            η_t_all[:,i] = cumtrapz(dat.t[what_pos],acf_η_tmp).*factor

            # Bulk viscosity
            acf_η_V_tmp = mean(acf_P[7:9])
            acf_η_V_all[:,i] = acf_η_V_tmp
            η_V_t_all[:,i] = cumtrapz(dat.t[what_pos],acf_η_V_tmp).*factor
        end
        what_pos1 = pos_step_start[1]:pos_step_end[1]

        # Cut viscosity due to large noise at long times
        if (mode == "single") & (mode_acf == "fft")
            cut_ratio = 0.9
            what_pos1 = what_pos1[1:round(Int64,length(what_pos1)*cut_ratio)]
            η_t_all = η_t_all[what_pos1,:]
            acf_η_all = acf_η_all[what_pos1,:]
            η_V_t_all = η_V_t_all[what_pos1,:]
            acf_η_V_all = acf_η_V_all[what_pos1,:]
        end

        if mode == "single"
            calc_error = 1
        elseif mode == "tdm"
            calc_error = 0
        end

        # Shear viscosity
        val, std, err = calc_average_GK(dat.step[what_pos1], η_t_all, info; do_err=calc_error, sym="η", name="viscosity", unit="Pa*s")
        plot_acf(dat.step[what_pos1],acf_η_all, info; sym="J_{η}^{(acf)}", name="viscosity", unit="Pa^2")
        η = single_dat(val, std, err)

        # Bulk viscosity
        val, std, err = calc_average_GK(dat.step[what_pos1], η_V_t_all, info; do_plt=0, do_err=0)
        η_V = single_dat(val, std, err)

        # Write data to file
        η_t = mean(η_t_all,dims=2)
        acf_η = mean(acf_η_all,dims=2)
        η_V_t = mean(η_V_t_all,dims=2)
        acf_η_V = mean(acf_η_V_all,dims=2)

        line1 = string("# Created by MD - Bulk Evaluation, Folder: ", info.folder)
        line2 = string("# t η(t) ACF_η η_V(t) ACF_η_V\n",
                       "# ps Pa*s Pa^2 Pa*s Pa^2")
        if (reduced_units) line2 = "# t* η*(t) ACF_η* η_V*(t) ACF_η_V*" end
        header = string(line1,"\n",line2)
        file = string(info.folder,"/viscosity.dat")
        fID = open(file,"w"); println(fID,header)
        writedlm(fID, hcat(dat.t[what_pos1],η_t,acf_η,η_V_t,acf_η_V), " ")
        close(fID)

        close("all")
    else
        η = single_dat(NaN,NaN,NaN)
        η_V = single_dat(NaN,NaN,NaN)
    end

    return η, η_V, T, ρ
end

## Thermal conducitvity --------------------------------------------------------
# Calculation of thermal conducitvity from heat current vector
function calc_thermalconductivity(info::info_struct, mode::String; mode_acf::String, CorrLength::Int64=0, SpanCorrFun::Int64=1, nEvery::Int64=1)
    # Loading Heat Flux File
    dat = load_heatflux(info,nEvery)
    if (nEvery > 1) @warn("Only every $(nEvery). step of heat flux vector used for ACF calculation.") end

    if (dat isa heat_dat)
        what = (dat.step .>= info.n_equ)

        # Unit conversion and GK factor
        metal2WmK = 2.566969967e-16
        factor = metal2WmK * mean(dat.V[what]) / (3 * kB * mean(dat.T[what])^2)
        if (reduced_units) factor = mean(dat.V[what]) / (3 * mean(dat.T[what])^2) end

        # Define start and end steps of individual correlation functions
        pos_step_start = []
        pos_step_end = []
        if mode == "single"
            for i_step in info.n_equ:SpanCorrFun:(maximum(dat.step) - CorrLength)
                append!(pos_step_start,findfirst(dat.step .>= i_step))
                append!(pos_step_end,findfirst(dat.step .>= i_step+CorrLength))
            end
        elseif mode == "tdm"
            append!(pos_step_start,findfirst(what))
            append!(pos_step_end,length(dat.step))
        end

        # Allocate arrays
        λ_t_all =   NaN .* ones(pos_step_end[1] .- pos_step_start[1]+1, length(pos_step_start))
        acf_λ_all = NaN .* ones(pos_step_end[1] .- pos_step_start[1]+1, length(pos_step_start))

        for i in 1:length(pos_step_start)
            what_pos = pos_step_start[i]:pos_step_end[i]

            # Calculate autocorrelation function
            J_mat = [dat.jx[what_pos]./dat.V[what_pos],
                     dat.jy[what_pos]./dat.V[what_pos],
                     dat.jz[what_pos]./dat.V[what_pos]]
            acf_J = pmap(x -> calc_acf(x, mode_acf),J_mat)

            acf_λ_tmp = (acf_J[1] .+ acf_J[2] .+ acf_J[3])
            acf_λ_all[:,i] = acf_λ_tmp
            λ_t_all[:,i] = cumtrapz(dat.t[what_pos],acf_λ_tmp).*factor
        end
        what_pos1 = pos_step_start[1]:pos_step_end[1]

        # Cut viscosity due to large noise at long times
        if (mode == "single") & (mode_acf == "fft")
            cut_ratio = 0.9
            what_pos1 = what_pos1[1:round(Int64,length(what_pos1)*cut_ratio)]
            λ_t_all = λ_t_all[what_pos1,:]
            acf_λ_all = acf_λ_all[what_pos1,:]
        end

        if mode == "single"
            calc_error = 1
        elseif mode == "tdm"
            calc_error = 0
        end

        # Averaging and error estimation
        val, std, err = calc_average_GK(dat.step[what_pos1], λ_t_all, info; do_err=calc_error, sym="λ", name="thermalconductivity", unit="W/(m*K)")
        plot_acf(dat.step[what_pos1],acf_λ_all, info; sym="J_{λ}^{(acf)}", name="thermalconductivity", unit="eV^2/(Å^4*ps)")
        λ = single_dat(val, std, err)

        # Write data to file
        λ_t = mean(λ_t_all,dims=2)
        acf_λ = mean(acf_λ_all,dims=2)
        line1 = string("# Created by MD - Bulk Evaluation, Folder: ", info.folder)
        line2 = string("# t λ(t) ACF_λ\n",
                       "# ps W/(m*K) eV^2/(Å^4*ps)")
        if (reduced_units) line2 = "# t* λ*(t) ACF_λ*" end
        header = string(line1,"\n",line2)
        file = string(info.folder,"/thermalconductivity.dat")
        fID = open(file,"w"); println(fID,header)
        writedlm(fID, hcat(dat.t[what_pos1],λ_t,acf_λ), " ")
        close(fID)

        close("all")
    else
        λ = single_dat(NaN,NaN,NaN)
    end

    return λ
end


## Subfunctions
# Function to integrate by trapezodial rule
@everywhere function cumtrapz(x::Array{Float64},y::Array{Float64})
    # Check input
    lx = length(x)
    ly = length(y)
    if (lx != ly) error("Length of input vectors not equal") end

    int = zeros(size(x))
    int[1] = 0
    for i=2:lx
        int[i] = int[i-1] + y[i-1] + (y[i]-y[i-1])/2
    end
    dx = sum(x[2:end] .- x[1:end-1])/(lx-1)
    int = int .* dx
    return int
end

# Function to calculate autocorrelation function by:
#   1. julia function "autocov"
#   2. FFT
@everywhere function calc_acf(in::Array{Float64,1},mode::String)
    L = length(in)
    if mode == "autocov"            # Mode: autocov
        out = autocov(in,Array(0:L-1),demean=false)
    elseif mode == "fft"            # Mode: FFT
        in0 = vcat(in,zeros(L))
        out = real(ifft(fft(in0) .* conj(fft(in0))))[1:L]
        out = out ./ (L:-1:1)
    end
    return out
end

# Function to calculate averge from GK integral
function calc_average_GK(steps, ave_t_all, info; do_plt=1, do_err=1, N_block=100, sym="", name="", unit="")
    ## Average value calculation
    # Running integral average
    ave_t = mean(ave_t_all,dims=2)[:]

    ## Fitting procedure following (2) Eq. 10
    # Function to fit following Ref. (3)
    @. fun_ave(x, p) = p[1] .* ( ((p[2].*p[3].*(1 .-exp(-x./p[3])) .+
                                  (1 .-p[1]).*p[4]).*(1 .-exp(-x./p[4]))) ./
                                  (p[2].*p[3].+(1 .-p[1]).*p[4]) )

    # Fitting
    p0 = [  [abs(ave_t[end]), 1.0, 1.0,  1.0],
            [abs(ave_t[end]), 1.0, 0.1,  10.0],
            [abs(ave_t[end]), 1.0, 1e-3, 1.0],
            [abs(ave_t[end]), 0.1, 0.1,  1.0],
            [0.1,             1.0, 0.1,  1.0] ]
    fit_ave = []
    converged = false
    k = 0
    while !(converged) && k < length(p0)
        k += 1
        fit_ave = curve_fit(fun_ave, steps, ave_t, p0[k])
        converged = fit_ave.converged
    end

    if fit_ave.converged
        val = fit_ave.param[1]
    else
        @warn("Fit to single data not converged!")
        val = NaN
    end

    ## Error bar calculation
    if do_err == 1
        # Get number of independent ACF's
        N = size(ave_t_all, 2)
        M_block = floor(Int64,N/N_block)

        M_block_min = 3
        if M_block < M_block_min
            std = NaN
            err = NaN
            @warn("M_block [$M_block] < M_block_min [$M_block_min]\nError calculation by block averaging not reasonable!")
        else
            # Loop all blocks
            vals = []
            for i = 1:N_block
                col_start = (i-1)*M_block + 1
                if i < N_block
                    what = col_start:(col_start+M_block-1)
                else
                    what = col_start:N
                end
                push!(vals,calc_average_GK(steps, ave_t_all[:,what], info; do_plt = 0, do_err = 0)[1])
            end

            # Exclude very high values
            vals = vals[abs.(vals) .< 4*val]

            # Calculation of standard deviation and standard error
            μ = mean(vals)
            std = sqrt(sum((vals.-μ).^2)/length(vals))
            err = std./sqrt(length(vals))
        end
    else
        std = NaN
        err = NaN
    end

    # Plotting
    if do_plt == 1
        # Reduce data if big dataset
        teve = ceil(Int64,size(ave_t_all,1)/10000)

        t = steps.*info.dt
        figure()
        plot(t[1:teve:end], ave_t[1:teve:end], "b", linewidth=1, label="running integral")
        plot(t,fun_ave(steps, fit_ave.param),"r", lw=0.75, label="fit to running integral")
        # fill_between(t,minimum(ave_t_all,dims=2)[:],maximum(ave_t_all,dims=2)[:],color="blue",alpha=0.2)
        plot([t[1],t[end]],[val,val],"k",label="result value", linewidth=0.5)
        if do_err == 1
            plot([t[1],t[end]], [val+std, val+std], "k--", linewidth=0.5, label="calculated standard deviation")
            plot([t[1],t[end]], [val-std, val-std], "k--", linewidth=0.5)
        end
        if !(reduced_units) xlabel("\$t\$ / ps") elseif reduced_units xlabel("\$t*\$") end
        if !(reduced_units) ylabel("\$$sym\$ / $unit") elseif reduced_units ylabel("\$$sym*\$") end
        legend(loc="lower right")
        tight_layout()
        savefig(string(info.folder,"/fig_$name(t).pdf"))

        close("all")
    end

    return val, std, err
end

# Function to calculate error bars from all GK integrals
function plot_acf(steps,acf_all, info; sym="", name="", unit="")
    # Reduce data if big dataset
    teve = ceil(Int64,size(acf_all,1)/10000)

    # Calculate average of acf
    t = steps.*info.dt
    ave_acf = mean(acf_all,dims=2)
    min_acf = minimum(acf_all,dims=2)[:]
    max_acf = maximum(acf_all,dims=2)[:]

    # Plots
    figure()
    plot(t[1:teve:end],ave_acf[1:teve:end], "b")
    fill_between(t[1:teve:end],min_acf[1:teve:end],max_acf[1:teve:end], color = "blue", alpha = 0.2)
    if !(reduced_units) xlabel("\$t\$ / ps") elseif reduced_units xlabel("\$t*\$") end
    if !(reduced_units) ylabel("\$$sym\$ / $unit") elseif reduced_units ylabel("\$$sym*\$") end
    tight_layout()
    savefig(string(info.folder,"/fig_$name-acf(t).pdf"))
    close("all")
end

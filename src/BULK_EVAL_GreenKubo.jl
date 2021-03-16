# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - GreenKubo
# Containing functions to apply Green-Kubo method to viscosity and thermal
# conducitvity
# ---
# created by Sebastian Schmitt, 13.03.2021
# ------------------------------------------------------------------------------
# References:
# (1) Humbert, M. T.; Zhang, Y.; Maginn, E. J. PyLAT: Python LAMMPS Analysis Tools. J. Chem. Inf. Model. 2019, 59 (4), 1301–1305. https://doi.org/10.1021/acs.jcim.9b00066.

## Settings --------------------------------------------------------------------
# Calculation of autocorrelation function (mainly in TDM method)
#   1 - autocov (julia function from StatsBase.jl -> slower, but "real" ACF)
#   2 - FFT (calculation by FFT algorithm -> very fast, but just estimation (??))
mode_acf_single = 1             # recommended: mode_acf_single = 1
mode_acf_tdm = 2                # mode_acf_tdm = 1 can be ver slowly

## Viscosity -------------------------------------------------------------------
# Calculation of viscosity from pressure tensor
# Mode: "single_run"
# -> calculation of ACF and η from divided parts and averaging/ error estimation
function calc_viscosities(info::info_struct, CorrLength::Int64, SpanCorrFun::Int64)
    # Loading Pressure File
    dat = load_pressure(info)

    if (dat isa pressure_dat)
        what = (dat.step .>= info.n_equ)

        # Unit conversion and GK factor
        metal2Pas = 1e-32
        factor = metal2Pas .* mean(dat.V[what]) ./ (kB .* mean(dat.T[what]))
        if (reduced_units) factor = mean(dat.V[what]) ./ mean(dat.T[what]) end

        # Define start and end steps of individual correlation functions
        pos_step_start = []
        pos_step_end = []
        for i_step in info.n_equ:SpanCorrFun:(maximum(dat.step) - CorrLength)
            append!(pos_step_start,findfirst(dat.step .>= i_step))
            append!(pos_step_end,findfirst(dat.step .>= i_step+CorrLength))
        end

        # Allocate arrays
        η_t_all = NaN .* ones(pos_step_end[1].-pos_step_start[1]+1,length(pos_step_start))
        η_V_t_all = NaN .* ones(pos_step_end[1].-pos_step_start[1]+1,length(pos_step_start))
        acf_η_all = NaN .* ones(pos_step_end[1].-pos_step_start[1]+1,length(pos_step_start))
        acf_η_V_all = NaN .* ones(pos_step_end[1].-pos_step_start[1]+1,length(pos_step_start))

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

            acf_P = pmap(x -> calc_acf(x,mode_acf_single),vcat(P_mat,P_diag))

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
        η_t = mean(η_t_all,dims=2)
        acf_η = mean(acf_η_all,dims=2)
        η_V_t = mean(η_V_t_all,dims=2)
        acf_η_V = mean(acf_η_V_all,dims=2)

        # Averaging and error estimation
        what_ave = dat.step[dat.step .<= CorrLength] .>= minimum([SpanCorrFun,round(CorrLength/3)])
        η = single_dat( mean(η_t[what_ave]),                                    # Mean value
                        block_average(η_t[what_ave])[1],                        # Block average
                        maximum(abs.(η_t[what_ave].-mean(η_t[what_ave]))))      # Max deviation
        η_V = single_dat(mean(η_V_t[what_ave]),                                 # Mean value
                        block_average(η_V_t[what_ave])[1],                      # Block average
                        maximum(abs.(η_V_t[what_ave].-mean(η_V_t[what_ave]))))  # Max deviation

        # Write data to file
        what_pos1 = pos_step_start[1]:pos_step_end[1]
        line1 = string("# Created by MD - Bulk Evaluation, Folder: ", info.folder)
        line2 = string("# t η(t) ACF_η η_V(t) ACF_η_V\n",
                       "# ps Pa*s Pa^2 Pa*s Pa^2")
        if (reduced_units) line2 = "# t* η*(t) ACF_η* η_V*(t) ACF_η_V*" end
        header = string(line1,"\n",line2)
        file = string(info.folder,"/viscosity.dat")
        fID = open(file,"w"); println(fID,header)
        writedlm(fID, hcat(dat.t[what_pos1],η_t,acf_η,η_V_t,acf_η_V), " ")
        close(fID)

        # Plots
        figure()
        plot(dat.t[what_pos1],acf_η, linewidth=0.5)
        if !(reduced_units) xlabel(L"t / ps") elseif reduced_units xlabel(L"t*") end
        if !(reduced_units) ylabel(L"|ACF| / Pa^2") elseif reduced_units ylabel(L"|ACF*|") end
        tight_layout()
        savefig(string(info.folder,"/fig_eta-acf(t).pdf"))

        figure()
        plot(dat.t[what_pos1], η_t, "b", linewidth=1)
        plot(dat.t[what_pos1][what_ave], ones(sum(what_ave)).*η.val, "k", linewidth=1, label="mean")
        plot(dat.t[what_pos1][what_ave], ones(sum(what_ave)).*(η.val + η.std), "k--", linewidth=0.5, label="standard deviation")
        plot(dat.t[what_pos1][what_ave], ones(sum(what_ave)).*(η.val - η.std), "k--", linewidth=0.5)
        plot(dat.t[what_pos1][what_ave], ones(sum(what_ave)).*maximum(η_t[what_ave]), "k:", linewidth=0.5, label="max/min")
        plot(dat.t[what_pos1][what_ave], ones(sum(what_ave)).*minimum(η_t[what_ave]), "k:", linewidth=0.5)
        if !(reduced_units) xlabel(L"t / ps") elseif reduced_units xlabel(L"t*") end
        if !(reduced_units) ylabel(L"η / Pa^2") elseif reduced_units ylabel(L"η*") end
        tight_layout()
        savefig(string(info.folder,"/fig_eta(t).pdf"))

        close("all")
    else
        η = single_dat(NaN,NaN,NaN)
        η_V = single_dat(NaN,NaN,NaN)
    end

    return η, η_V
end

# mode: "tdm"
# -> calculation of ACF and η from total trajectory for TDM method
function calc_viscosities(info::info_struct)
    # Loading Pressure File
    dat = load_pressure(info)

    if (dat isa pressure_dat)
        what = (dat.step .>= info.n_equ)
        t = dat.t[what]

        # Calculate autocorrelation function
        P_mat = [dat.pxy[what],
                 dat.pxz[what],
                 dat.pyz[what],
                 dat.pxx[what] .- dat.pyy[what],
                 dat.pyy[what] .- dat.pzz[what],
                 dat.pxx[what] .- dat.pzz[what]]
        acf_P = pmap(x -> calc_acf(x,mode_acf_tdm),P_mat)

        # Unit conversion and GK factor
        metal2Pas = 1e-32
        factor = metal2Pas .* dat.V[what] ./ (kB .* dat.T[what])
        if (reduced_units) factor = dat.V[what] ./ (dat.T[what]) end

        # Integration of acf
        # -> method adopted from PyLat code (1)
        acf_η = (acf_P[1] .+ acf_P[2] .+ acf_P[3])./6 .+
                (acf_P[4] .+ acf_P[5] .+ acf_P[6])./24
        η_t = cumtrapz(t,acf_η).*factor

        # Bulk viscosity
        P_diag = [dat.pxx[what], dat.pyy[what], dat.pzz[what]]
        acf_η_V = mean(pmap(x -> calc_acf(x,mode_acf_tdm),P_diag))
        η_V_t = cumtrapz(t,acf_η_V).*factor

        # Write data to file
        line1 = string("# Created by MD - Bulk Evaluation, Folder: ", info.folder)
        line2 = string("# t η(t) ACF_η η_V(t) ACF_η_V\n",
                       "# ps Pa*s Pa^2 Pa*s Pa^2")
        if (reduced_units) line2 = "# t* η*(t) ACF_η* η_V*(t) ACF_η_V*" end
        header = string(line1,"\n",line2)
        file = string(info.folder,"/viscosity.dat")
        fID = open(file,"w"); println(fID,header)
        writedlm(fID, hcat(dat.t[what],η_t,acf_η,η_V_t,acf_η_V), " ")
        close(fID)

        # Plots
        teve = 100        # Just take data every 'teve' timesteps
        figure()
        plot(t[2:teve:end],abs.(acf_η[2:teve:end]), linewidth=0.1)
        loglog()
        if !(reduced_units) xlabel(L"t / ps") elseif reduced_units xlabel(L"t*") end
        if !(reduced_units) ylabel(L"|ACF| / Pa^2") elseif reduced_units ylabel(L"|ACF*|") end
        tight_layout()
        savefig(string(info.folder,"/fig_eta-acf(t).pdf"))

        figure()
        plot(t[1:teve:end],η_t[1:teve:end], linewidth=0.1)
        if !(reduced_units) xlabel(L"t / ps") elseif reduced_units xlabel(L"t*") end
        if !(reduced_units) ylabel(L"η / Pa^2") elseif reduced_units ylabel(L"η*") end
        tight_layout()
        savefig(string(info.folder,"/fig_eta(t).pdf"))
        close()

        # Save Results
        η = single_dat(η_t[end], NaN, NaN)
        η_V = single_dat(η_V_t[end], NaN, NaN)
    else
        η = single_dat(NaN,NaN,NaN)
        η_V = single_dat(NaN,NaN,NaN)
    end

    return η, η_V
end

## Thermal conducitvity --------------------------------------------------------
# Calculation of thermal conducitvity from heat current vector

# Mode: "single_run"
# -> calculation of ACF and λ from divided parts and averaging/ error estimation
function calc_thermalconductivity(info::info_struct, CorrLength::Int64, SpanCorrFun::Int64)
    # Loading Heat Flux File
    dat = load_heatflux(info)

    if (dat isa heat_dat)
        what = (dat.step .>= info.n_equ)

        # Unit conversion and GK factor
        metal2WmK = 2.566969967e-16
        factor = metal2WmK * mean(dat.V[what]) / (3 * kB * mean(dat.T[what])^2)
        if (reduced_units) factor = mean(dat.V[what]) / (3 * mean(dat.T[what])^2) end

        # Define start and end steps of individual correlation functions
        pos_step_start = []
        pos_step_end = []
        for i_step in info.n_equ:SpanCorrFun:(maximum(dat.step) - CorrLength)
            append!(pos_step_start,findfirst(dat.step .>= i_step))
            append!(pos_step_end,findfirst(dat.step .>= i_step+CorrLength))
        end

        # Allocate arrays
        λ_t_all = NaN .* ones(pos_step_end[1].-pos_step_start[1]+1,length(pos_step_start))
        acf_λ_all = NaN .* ones(pos_step_end[1].-pos_step_start[1]+1,length(pos_step_start))

        for i in 1:length(pos_step_start)
            what_pos = pos_step_start[i]:pos_step_end[i]

            # Calculate autocorrelation function
            J_mat = [dat.jx[what_pos]./dat.V[what_pos],
                     dat.jy[what_pos]./dat.V[what_pos],
                     dat.jz[what_pos]./dat.V[what_pos]]
            acf_J = pmap(x -> calc_acf(x,mode_acf_single),J_mat)

            acf_λ_tmp = (acf_J[1] .+ acf_J[2] .+ acf_J[3])
            acf_λ_all[:,i] = acf_λ_tmp
            λ_t_all[:,i] = cumtrapz(dat.t[what_pos],acf_λ_tmp).*factor
        end
        λ_t = mean(λ_t_all,dims=2)
        acf_λ = mean(acf_λ_all,dims=2)

        # Averaging and error estimation
        what_ave = dat.step[dat.step .<= CorrLength] .>= minimum([SpanCorrFun,round(CorrLength/3)])
        λ = single_dat( mean(λ_t[what_ave]),                                    # Mean value
                        block_average(λ_t[what_ave])[1],                        # Block average
                        maximum(abs.(λ_t[what_ave].-mean(λ_t[what_ave]))))      # Max deviation

        # Write data to file
        what_pos1 = pos_step_start[1]:pos_step_end[1]
        line1 = string("# Created by MD - Bulk Evaluation, Folder: ", info.folder)
        line2 = string("# t λ(t) ACF_λ\n",
                       "# ps W/(m*K) eV^2/(Å^4*ps)")
        if (reduced_units) line2 = "# t* λ*(t) ACF_λ*" end
        header = string(line1,"\n",line2)
        file = string(info.folder,"/thermalconductivity.dat")
        fID = open(file,"w"); println(fID,header)
        writedlm(fID, hcat(dat.t[what_pos1],λ_t,acf_λ), " ")
        close(fID)

        # Plots
        figure()
        plot(dat.t[what_pos1],acf_λ, linewidth=0.5)
        if !(reduced_units) xlabel(L"t / ps") elseif reduced_units xlabel(L"t*") end
        if !(reduced_units) ylabel(L"|ACF| / eV^2/(Å^4*ps)") elseif reduced_units ylabel(L"|ACF*|") end
        tight_layout()
        savefig(string(info.folder,"/fig_lambda-acf(t).pdf"))

        figure()
        plot(dat.t[what_pos1], λ_t, "b", linewidth=1)
        plot(dat.t[what_pos1][what_ave], ones(sum(what_ave)).*λ.val, "k", linewidth=1, label="mean")
        plot(dat.t[what_pos1][what_ave], ones(sum(what_ave)).*(λ.val + λ.std), "k--", linewidth=0.5, label="standard deviation")
        plot(dat.t[what_pos1][what_ave], ones(sum(what_ave)).*(λ.val - λ.std), "k--", linewidth=0.5)
        plot(dat.t[what_pos1][what_ave], ones(sum(what_ave)).*maximum(λ_t[what_ave]), "k:", linewidth=0.5, label="max/min")
        plot(dat.t[what_pos1][what_ave], ones(sum(what_ave)).*minimum(λ_t[what_ave]), "k:", linewidth=0.5)
        if !(reduced_units) xlabel(L"t / ps") elseif reduced_units xlabel(L"t*") end
        if !(reduced_units) ylabel(L"λ / (W/(m*K)") elseif reduced_units ylabel(L"λ*") end
        tight_layout()
        savefig(string(info.folder,"/fig_lambda(t).pdf"))

        close("all")
    else
        λ = single_dat(NaN,NaN,NaN)
    end

    return λ
end

# Mode: "tdm"
# -> calculation of ACF and λ from total trajectory for TDM method
function calc_thermalconductivity(info::info_struct)
    # Loading Heat Flux File
    dat = load_heatflux(info)

    if (dat isa heat_dat)
        what = (dat.step .>= info.n_equ)
        t = dat.t[what]

        # Calculate autocorrelation function
        J_mat = [dat.jx[what]./dat.V[what],
                 dat.jy[what]./dat.V[what],
                 dat.jz[what]./dat.V[what]]
        acf_J = pmap(x -> calc_acf(x,mode_acf_tdm),J_mat)

        # Unit conversion and GK factor
        metal2WmK = 2.566969967e-16
        factor = metal2WmK .* dat.V[what] ./ (3 .* kB .* dat.T[what].^2)
        if (reduced_units) factor = dat.V[what] ./ (3 .* dat.T[what].^2) end

        # Integration of acf
        acf_λ = (acf_J[1] .+ acf_J[2] .+ acf_J[3])
        λ_t = cumtrapz(t,acf_λ).*factor

        # Write data to file
        line1 = string("# Created by MD - Bulk Evaluation, Folder: ", info.folder)
        line2 = string("# t λ(t) ACF_λ\n",
                       "# ps W/(m*K) eV^2/(Å^4*ps)")
        if (reduced_units) line2 = "# t* λ*(t) ACF_λ*" end
        header = string(line1,"\n",line2)
        file = string(info.folder,"/thermalconductivity.dat")
        fID = open(file,"w"); println(fID,header)
        writedlm(fID, hcat(dat.t[what],λ_t,acf_λ), " ")
        close(fID)

        # Plots
        teve = 100        # Just take data every 'teve' timesteps
        figure()
        tight_layout()
        abs_λ = abs.(acf_λ[2:teve:end])
        plot(t[2:teve:end],abs_λ, linewidth=0.1)
        loglog()
        if !(reduced_units) xlabel(L"t / ps") elseif reduced_units xlabel(L"t*") end
        if !(reduced_units) ylabel(L"|ACF| / eV^2/(Å^4*ps)") elseif reduced_units ylabel(L"|ACF*|") end
        tight_layout()
        savefig(string(info.folder,"/fig_lambda_acf(t).pdf"))
        close()

        figure()
        plot(t[1:teve:end],λ_t[1:teve:end], linewidth=0.1)
        if !(reduced_units) xlabel(L"t / ps") elseif reduced_units xlabel(L"t*") end
        if !(reduced_units) ylabel(L"λ / (W/(m*K)") elseif reduced_units ylabel(L"λ*") end
        tight_layout()
        savefig(string(info.folder,"/fig_lambda(t).pdf"))
        close()

        λ = single_dat(λ_t[end],NaN,NaN)
    else
        λ = single_dat(NaN,NaN,NaN)
    end

    return λ
end

## General functions
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
@everywhere function calc_acf(in::Array{Float64,1},mode::Int64)
    L = length(in)
    if mode == 1            # Mode: autocov
        out = autocov(in,Array(0:L-1),demean=false)
    elseif mode == 2        # Mode: FFT
        in0 = vcat(in,zeros(L))
        out = real(ifft(fft(in0) .* conj(fft(in0))))[1:L]
        out = out ./ (L:-1:1)
    end
    return out
end

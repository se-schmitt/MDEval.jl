## NEMDShear.jl
# ------------------------------------------------------------------------------
# Evaluation Software for NEMD shear simulations (with sllod algorithm)
# Containing functions to evaluate output files of NEMD shear simulations to
# calculate viscosity
# ---
# created by Sebastian Schmitt, 16.07.2021

function EvalNEMDShear(subfolder,inpar)
    moltype, dt, natoms, molmass = load_info(subfolder)

    # Initialization of info structure
    info = info_struct( subfolder,          # info.folder
                        inpar.ensemble,     # info.ensemble
                        inpar.n_equ,        # info.n_equ
                        moltype,            # info.moltype
                        dt,                 # info.dt
                        natoms,             # info.natoms
                        molmass,            # info.molmass
                        inpar.n_blocks,     # info.n_blocks
                        inpar.N_bin,        # info.N_bin
                        inpar.r_cut)        # info.r_cut

    # Average Thermodynamic Properties
    T, p, ρ, Etot, Ekin, Epot, c, pyz = ave_thermo(info; is_nemd="shear")

    # Load profile data
    list = readdir(info.folder)
    ts_add = 0
    prof = []
    no_chunks = 0
    for f in list[startswith.(list,"vy_profile")]
        prof, no_chunks, n_steps = read_profile1D("$(info.folder)/$f",prof,ts_add)
        ts_add = prof.timestep[end]
    end
    if info.n_equ > 0
        what_eval = prof.timestep .≥ info.n_equ
        prof.timestep = prof.timestep[what_eval]
        for f in fieldnames(typeof(prof))
            if !(f in [:timestep])
                setfield!(prof,f,getfield(prof,f)[what_eval,:])
            end
        end
        prof.timestep = prof.timestep .- prof.timestep[1]
    end

    # Load thermo data
    dat = load_thermo(info, is_nemd="shear")
    if info.n_equ > 0
        what_eval = dat.step .≥ info.n_equ
        for f in fieldnames(typeof(dat))
            if !(f in [:method]) && !isempty(getfield(dat,f))
                setfield!(dat,f,getfield(dat,f)[what_eval])
            end
        end
        dat.step = dat.step .- dat.step[1]
        dat.t = dat.t .- dat.t[1]
    end
    if (reduced_units)      factor_p = 1; factor_t = 1
    elseif !(reduced_units) factor_p = 1e5; factor_t = 10^12 end

    Lx, Ly, Lz = load_info(info.folder; is_nemd="shear")

    # Calculate velocity gradient
    xbins = mean(prof.x,dims=1)[:] .* Lz
    vybins = mean(prof.vy,dims=1)[:]
    vybins_std = std(prof.vy,dims=1)[:]

    if dat.method == "sllod"
        # SLLOD
        # Fit velocity profile
        fit_v = ([ones(length(xbins[2:end-1]),1) xbins[2:end-1]] \ vybins[2:end-1])
        dvdx = fit_v[2]
        yinterc = fit_v[1]
        s_rate = factor_t*dvdx

        # Calculation of viscosity and its uncertainties
        η_vec = -dat.pyz.*factor_p/s_rate
        η_ave = mean(η_vec)
        η_std, η_err = block_average(η_vec; N_blocks=100)

    elseif dat.method == "rnemd"
        # RNEMD
        N_blocks = 5

        # Fit velocity profiles
        what_grad1 = 2:round(Int64,no_chunks/2)
        what_grad2 = round(Int64,no_chunks/2+2):no_chunks
        dvdx_all = Float64[]
        yinterc_all = Float64[]
        for i in 1:N_blocks
            what_i = round(Int64,size(prof.x,1)/N_blocks*(i-1)+1):round(Int64,size(prof.x,1)/N_blocks*i)
            (α1_i,dvdx1_i) = [ones(length(what_grad1)) mean(prof.x[what_i,:].*Lz,dims=1)[what_grad1]] \ mean(prof.vy[what_i,:],dims=1)[what_grad1]
            (α2_i,dvdx2_i) = [ones(length(what_grad2)) mean(prof.x[what_i,:].*Lz,dims=1)[what_grad2]] \ mean(prof.vy[what_i,:],dims=1)[what_grad2]
            append!(dvdx_all,abs.([dvdx1_i,dvdx2_i]))
            append!(yinterc_all,[α1_i,α2_i])
        end
        s_rate = factor_t * mean(dvdx_all)
        dvdx1 = mean(dvdx_all[1:2:end])
        dvdx2 = mean(dvdx_all[2:2:end])
        yinterc1 = mean(yinterc_all[1:2:end])
        yinterc2 = mean(yinterc_all[2:2:end])
        s_rate_std = factor_t * std(dvdx_all)

        # Calculate slope to accumulated stress
        pyz_all = Float64[]
        for i in 1:N_blocks
            what_i = round(Int64,length(dat.t)/N_blocks*(i-1)+1):round(Int64,length(dat.t)/N_blocks*i)
            (α_i,pyz_i) = [ones(length(what_i)) dat.t[what_i]] \ dat.pyz[what_i]
            push!(pyz_all,pyz_i)
        end
        pyz = mean(pyz_all)
        pyz_std = std(pyz_all)
        
        # Calculation of viscosity and its uncertainties
        η_ave = -pyz*0.5.*factor_p/s_rate
        η_std = -pyz*0.5.*factor_p/s_rate^2 * s_rate_std + pyz_std*0.5.*factor_p/s_rate
        η_err = η_std * quantile(TDist(length(xbins)-2),0.975)
        
    end

    # Save viscosity
    η = single_dat(η_ave, η_std, η_err)

    # Output results
    res = results_struct_nemd(T, p, ρ, [1.0], Etot, Ekin, Epot, η, s_rate, [])
    OutputResultNEMD(res, info.folder)

    # Figures ------------------------------------------------------------------
    # Figure of velocity profile
    figure()
    title("Mean velocity profile of simulation box")
    if reduced_units    xlabel(L"z^*");                 ylabel(L"v_y^*")
    else                xlabel(string(L"z"," / Å"));    ylabel(string(L"v_y"," / Å/ps"))    end
    errorbar(xbins, vybins, yerr=vybins_std, fmt="bo", capsize=2, mfc="none")
    if dat.method == "sllod"
        plot(xbins[2:end-1], xbins[2:end-1].*dvdx.+yinterc, color="red", label="Fit")
    elseif dat.method == "rnemd"
        plot(xbins[what_grad1], xbins[what_grad1].*dvdx1.+yinterc1, color="red", label="Fit gradient 1")
        plot(xbins[what_grad2], xbins[what_grad2].*dvdx2.+yinterc2, color="red", label="Fit gradient 2")
    end
    xlim([0,Lz])
    legend(loc="upper left")
    tight_layout()
    savefig("$(info.folder)/fig_vy_profile.pdf")

    # Dependency of the velocity gradient on time
    figure()
    w_sz = 20
    if dat.method == "sllod"
        plot(prof.timestep[Int64(w_sz/2):end-Int64(w_sz/2)].*info.dt, rollmean(mean(prof.vy[:,1],dims=2)[:],w_sz), "-b", label="\$0.0 \\leq z/L_{\\rm z} \\leq $(round(xbins[1],digits=4))\$")
        plot(prof.timestep[Int64(w_sz/2):end-Int64(w_sz/2)].*info.dt, rollmean(mean(prof.vy[:,end],dims=2)[:],w_sz), "-r", label="\$$(round(xbins[end],digits=4)) \\leq z/L_{\\rm z} \\leq 1.0\$")
    elseif dat.method == "rnemd"
        plot(prof.timestep[Int64(w_sz/2):end-Int64(w_sz/2)].*info.dt, rollmean(prof.vy[:,1][:],w_sz), "-b", label="bin 1 (\$v = v_{\\rm min}\$)")
        plot(prof.timestep[Int64(w_sz/2):end-Int64(w_sz/2)].*info.dt, rollmean(prof.vy[:,round(Int64,no_chunks/2+1)][:],w_sz), "-b", label="bin $(round(Int64,no_chunks/2+1)) (\$v = v_{\\rm max}\$)")
    end
    if reduced_units    xlabel("\$t^*\$");              ylabel("\$v_{\\rm y}^*\$")
    else                xlabel("\$t~/~{\\rm ps}\$");    ylabel("\$v_{\\rm y}~/~{\\rm Å~ps^{-1}}\$")    end
    legend()
    xlim([0,prof.timestep[end]*info.dt])
    tight_layout()
    savefig("$(info.folder)/fig_time_dependence_vy.pdf")

    close("all")
end

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
    for f in list[startswith.(list,"vy_profile")]
        prof, no_chunks, n_steps = read_profile1D("$(info.folder)/$f",prof,ts_add)
        ts_add = prof.timestep[end]
    end

    # Load thero data
    dat = load_thermo(info, is_nemd="shear")
    if (reduced_units)      factor_p = 1; factor_t = 1
    elseif !(reduced_units) factor_p = 1e5; factor_t = 10^12 end

    # Calculate velocity gradient
    xbins = mean(prof.x,dims=1)[:]
    vybins = mean(prof.vy,dims=1)[:]
    vybins_std = std(prof.vy,dims=1)[:]
    dvdx = ([ones(length(xbins[2:end-1]),1) xbins[2:end-1]] \ vybins[2:end-1])[2]
    yinterc = ([ones(length(xbins[2:end-1]),1) xbins[2:end-1]] \ vybins[2:end-1])[1]

    Lx, Ly, Lz = load_info(info.folder; is_nemd="shear")
    s_rate = factor_t*dvdx/Lz

    # Calculation of viscosity and its uncertainties
    η_vec = -dat.pyz.*factor_p/s_rate
    η_ave = mean(η_vec)
    η_std, η_err = block_average(η_vec; M_block=100)

    # Save viscosity
    η = single_dat(η_ave, η_std, η_err)

    # Output results
    res = results_struct_nemd(T, p, ρ, [1.0], Etot, Ekin, Epot, η, s_rate, [])
    OutputResultNEMD(res, info.folder)

    # Figures ------------------------------------------------------------------
    # Figure of velocity profile
    figure()
    title("Mean velocity profile of simulation box")
    if reduced_units    xlabel(L"z^* / L_z^*");                 ylabel(L"v_y^*")
    else                xlabel(string(L"z / L_z"," / Å"));    ylabel(string(L"v_y"," / Å/ps"))    end
    errorbar(xbins, vybins, yerr=vybins_std, fmt="bo-", capsize=2)
    plot(xbins[2:end-1], xbins[2:end-1].*dvdx.+yinterc, color="red", label="Fit")
    legend(loc="upper left")
    tight_layout()
    savefig("$(info.folder)/fit_vy_profile.pdf")
end

## EvalNEMDHeat.jl
# ------------------------------------------------------------------------------
# Evaluation Software for NEMD heat simulations
# Functions to analyse temperature profile and calculate thermal conductivity
# ---
# created by Jens Wagner, 04.01.2021
# ------------------------------------------------------------------------------

# EvalNEMDHeat
function EvalNEMDHeat(subfolder, inpar)

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

    # Average thermodynamic properties
    T, p, ρ, Etot, Ekin, Epot, c, _ = ave_thermo(info; is_nemd="heat")

    # Evaluate temperature profile
    list = readdir(info.folder)
    prof = []
    ts_add = 0
    for f in list[startswith.(list,"temp_profile")]
        prof, nbins, nts = read_profile1D("$(info.folder)/$f",prof,ts_add)
        ts_add = prof.timestep[end]
    end

    nbins = size(prof.x,2)
    xbins = dropdims(mean(prof.x,dims=1),dims=1)
    xbins = xbins .- xbins[1] .+ (xbins[2] - xbins[1])/2
    Tbins = dropdims(mean(prof.T,dims=1),dims=1)
    Tbins_std = dropdims(std(prof.T,dims=1),dims=1)

    # Read info file
    Lx,Ly,Lz = load_info(info.folder; is_nemd="heat")
    A = Ly*Lz

    # Evaluate heat profile
    dat = load_thermo(info, is_nemd="heat")

    dQhot = abs.(dat.Qhot .- dat.Qhot[1])
    dQcold = abs.(dat.Qcold .- dat.Qcold[1])

    if (minimum(xbins) > 0.0) && (maximum(xbins) < 1.0)
        xbins = xbins.*Lx
    end

    buffer = 0.02*Lx
    L_thermo = inpar.k_L_thermo*Lx
    what_zone1 = (xbins .> (L_thermo/2 + buffer))        .& (xbins .< (Lx/2 - L_thermo/2 - buffer))
    what_zone2 = (xbins .> (Lx/2 + L_thermo/2 + buffer)) .& (xbins .< (Lx - L_thermo/2 - buffer))

    # Calculate thermal conductivity -------------------------------------------
    λval, dQ_fit, fit_zone1, fit_zone2 = calc_lambda_NEMD(dat.t, dQhot, dQcold, xbins, Tbins, what_zone1, what_zone2, A)

    # Block averaging ----------------------------------------------------------
    λvec = Array{Float64,1}(undef,info.n_blocks)
    sz_dat = floor(Int64,length(dat.t)/info.n_blocks)
    sz_prof = floor(Int64,size(prof.T,1)/info.n_blocks)
    for i = 1:info.n_blocks
        what_i_dat = ((i-1)*sz_dat+1):(i*sz_dat)
        what_i_prof = ((i-1)*sz_prof+1):(i*sz_prof)
        λvec[i], _ , _ , _ = calc_lambda_NEMD(dat.t[what_i_dat], dQhot[what_i_dat], dQcold[what_i_dat], xbins, dropdims(mean(prof.T[what_i_prof,:],dims=1),dims=1), what_zone1, what_zone2, A)
    end
    λ = single_dat(λval,std(λvec),std(λvec)/sqrt(info.n_blocks))

    res = results_struct_nemd(T, p, ρ, [1.0], Etot, Ekin, Epot, [], [], λ)
    OutputResultNEMD(res,info.folder)

    # Figures ------------------------------------------------------------------
    # Figure of temperature profile
    figure()
    title("Temperature profile of simulation box")
    if reduced_units    xlabel(L"x^*");                 ylabel(L"T^*")
    else                xlabel(string(L"x"," / Å"));    ylabel(string(L"T"," / K"))    end
    what_not_zero = Tbins .> 0.0
    errorbar(xbins[what_not_zero], Tbins[what_not_zero], yerr=Tbins_std[what_not_zero], fmt="bo-", capsize=2, zorder=1)
    plot(xbins[what_zone1],fit_zone1[1].+xbins[what_zone1].*fit_zone1[2],"k-", zorder=2)
    plot(xbins[what_zone2],fit_zone2[1].+xbins[what_zone2].*fit_zone2[2],"k-", zorder=3)
    tight_layout()
    savefig(string(info.folder, "/T_profile.pdf"))

    # Figure of heats
    figure()
    title("Heat profile of simulation box")
    if reduced_units    xlabel(L"t*");                  ylabel(L"Q^*")
    else                xlabel(string(L"t"," / ps"));   ylabel(string(L"Q"," / eV"))    end
    plot(dat.t, dQhot, label="dQ_hot")
    plot(dat.t, dQcold, label="dQ_cold")
    plot(dat.t, dat.t.*dQ_fit, label="Fit")
    legend(loc="upper left");   tight_layout()
    savefig(string(info.folder, "/Q_profile.pdf"))
end

## Subfunctions ----------------------------------------------------------------
function calc_lambda_NEMD(t, dQhot, dQcold, xbins, Tbins, what_zone1, what_zone2, A)
    # Linear fit of heat
    dQ_fit = (repeat([ones(length(t)) t],2,1) \ [dQhot;dQcold])[2]

    # Linear fit to temperature profile
    fit_zone1 = [ones(sum(what_zone1)) xbins[what_zone1]] \ Tbins[what_zone1]
    dTdx1 = abs(fit_zone1[2])
    fit_zone2 = [ones(sum(what_zone2)) xbins[what_zone2]] \ Tbins[what_zone2]
    dTdx2 = abs(fit_zone2[2])
    dTdx = (dTdx1 + dTdx2)/2

    λ = dQ_fit / A / dTdx / 2

    return λ, dQ_fit, fit_zone1, fit_zone2
end

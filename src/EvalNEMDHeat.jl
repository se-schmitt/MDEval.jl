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

    # Evaluate heat profile
    dat = load_thermo(info, is_nemd="heat")

    dQhot = abs.(dat.Qhot .- dat.Qhot[1])
    dQcold = abs.(dat.Qcold .- dat.Qcold[1])

    # Linear regression to a linear funtion
    dQ_fit = (repeat([ones(length(dat.t)) dat.t],2,1) \ [dQhot;dQcold])[2]

    # Calculate thermal conductivity -------------------------------------------
    k_L_thermo = 0.25

    # Read info file
    Lx,Ly,Lz = add_info_NEMDheat(info)
    A = Ly*Lz

    if (minimum(xbins) > 0.0) && (maximum(xbins) < 1.0)
        xbins = xbins.*Lx
    end

    buffer = 0.05*Lx
    L_thermo = k_L_thermo*Lx
    what_zone1 = (xbins .> (L_thermo/2 + buffer))        .& (xbins .< (Lx/2 - L_thermo/2 - buffer))
    what_zone2 = (xbins .> (Lx/2 + L_thermo/2 + buffer)) .& (xbins .< (Lx - L_thermo/2 - buffer))

    fit_zone1 = [ones(sum(what_zone1)) xbins[what_zone1]] \ Tbins[what_zone1]
    dTdx1 = abs(fit_zone1[2])
    fit_zone2 = [ones(sum(what_zone2)) xbins[what_zone2]] \ Tbins[what_zone2]
    dTdx2 = abs(fit_zone2[2])
    dTdx = (dTdx1 + dTdx2)/2

    λ = dQ_fit / A / dTdx

    res = results_struct_nemd(T, p, ρ, [1.0], Etot, Ekin, Epot, [], [], [], [], λ)
    OutputResultNEMD(res,info.folder)
    @exfiltrate

    # Figures ------------------------------------------------------------------
    # Figure of temperature profile
    figure()
    title("Mean temperature profile of simulation box")
    xlabel("Bin")
    if reduced_units        ylabel(L"T^*")
    else                    ylabel(string(L"T"," / K"))    end
    errorbar(xbins, Tbins, yerr=Tbins_std, fmt="bo-", capsize=2)
    plot(xbins[what_zone1],fit_zone1[1].+xbins[what_zone1].*fit_zone1[2],"k-")
    plot(xbins[what_zone2],fit_zone2[1].+xbins[what_zone2].*fit_zone2[2],"k-")
    savefig(string(info.folder, "/T_profile.pdf"))

    # Figure of heats
    figure()
    title("Heat profile of simulation box")
    if reduced_units    xlabel(L"t*");                  ylabel(L"Q^*")
    else                xlabel(string(L"t"," / ps"));   ylabel(string(L"Q"," / eV"))    end
    plot(dat.t, dQhot, label="dQ_hot")
    plot(dat.t, dQcold, label="dQ_cold")
    plot(dat.t, dat.t.*dQ_fit, label="Fit")
    legend()
    savefig(string(info.folder, "/Q_profile.pdf"))
end

## Subfunctions ----------------------------------------------------------------
# Load additional info
function add_info_NEMDheat(info)

    # Get filename
    list = readdir(info.folder)
    fname = list[occursin.("info.", list)][end]
    file = string(info.folder,"/", fname)

    # Read file info.dat
    if isfile(file)
        fID = open(file,"r");   lines = readlines(fID);     close(fID)

        # Extract information
        pos5 = findfirst(": ",lines[5])
        Lx = parse(Float64,lines[5][pos5[end]+1:end])
        pos6 = findfirst(": ",lines[6])
        Ly = parse(Float64,lines[6][pos6[end]+1:end])
        pos7 = findfirst(": ",lines[7])
        Lz = parse(Float64,lines[7][pos7[end]+1:end])
        pos8 = findfirst(": ",lines[8])
        nbins = parse(Float64,lines[8][pos8[end]+1:end])

    else error("File \"",file,"\" is empty") end

    return Lx, Ly, Lz
end

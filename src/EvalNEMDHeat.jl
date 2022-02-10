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
    T, p, ρ, Etot, Ekin, Epot, c, _ = ave_thermo(info, is_nemd="heat")
    
    # Evaluate temperature profile
    nbins, dims = add_info_NEMDheat(info)
    bstps, bins = load_profile(info, nbins)
    
    mbins = []
    for j = 1:length(bins)
        append!(mbins, mean(bins[j]))
    end
    
    plot_profile(info, mbins)
    plot_bins(info, bstps, bins)

    # Evaluate heat profile
    dat = load_thermo(info, is_nemd="heat")

    dQhot = broadcast(abs, dat.Qhot .- dat.Qhot[1])
    dQcold = dat.Qcold .- dat.Qcold[1]
    dQ = []
    for j = 1:length(dQhot)
        append!(dQ, (dQhot[j] + dQcold[j]) /2)
    end
    _, fit_dQ = linear_fit(dat.t, dQ)

    plot_Q(info, dat.t, dQhot, dQcold, dQ, fit_dQ)

    # Calculate thermal conductivity
    λ = calc_λ(nbins,mbins,dims,fit_dQ)

    res = results_struct_nemd_heat(T, p, ρ, c, Etot, Ekin, Epot, λ)
    OutputResultNEMDHeat(res,info.folder)
end

# Subfunctions
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
        xlen = parse(Float64,lines[5][pos5[end]+1:end])
        pos6 = findfirst(": ",lines[6])
        ylen = parse(Float64,lines[6][pos6[end]+1:end])
        pos7 = findfirst(": ",lines[7])
        zlen = parse(Float64,lines[7][pos7[end]+1:end])
        pos8 = findfirst(": ",lines[8])
        nbins = parse(Float64,lines[8][pos8[end]+1:end])
        
    else error("File \"",file,"\" is empty") end
    
    dims = [xlen, ylen, zlen]
    return nbins, dims
end

# Load temperature profile
function load_profile(info, nbins)

    fname = "$(info.folder)/temp_profile.dat"
    fID = open(fname,"r");    lines = readlines(fID);    close(fID)

    bins = [[] for i in 1:nbins]
    bstps = []
    for i = 3:length(lines)
        for j = 1:trunc(Int, nbins)
            if startswith(lines[i],"  $(j) ")
                if j == 1
                    append!(bstps, parse(Float64, split(lines[i-1]," ")[1]))
                end    
                append!(bins[j], parse(Float64, split(lines[i]," ")[6]))
            end
        end
    end

    return bstps, bins
end

# Plot temperature profile of simulation box
function plot_profile(info, mbins)
    
    figure()
    title("Mean temperature profile of simulation box")
    
    xlabel("Bin")
    if reduced_units
        ylabel(L"T^*")
    else
        ylabel(string(L"T"," / K"))
    end

    plot(1:length(mbins), mbins, "bo-")
    savefig(string(info.folder, "/T_profile.pdf"))

end

# Plot temperature profile of bins
function plot_bins(info, bstps, bins)

    figure()
    title("Temperature profile of bins")

    if reduced_units
        xlabel(L"t*")
        ylabel(L"T^*")
    else
        xlabel(string(L"t"," / ps"))
        ylabel(string(L"T"," / K"))
    end

    for j = 1:length(bins)
        plot(bstps.*info.dt, bins[j], label="bin$(j)")
    end
    legend()
    savefig(string(info.folder, "/T_bins.pdf"))

end

# Plot heat profile
function plot_Q(info, t, dQhot, dQcold, dQ, fit_dQ)

    figure()
    title("Heat profile of simulation box")

    if reduced_units
        xlabel(L"t*")
        ylabel(L"Q^*")
    else
        xlabel(string(L"t"," / ps"))
        ylabel(string(L"Q"," / eV"))
    end

    plot(t, dQhot, label="dQ_hot")
    plot(t, dQcold, label="dQ_cold")
    plot(t, dQ, label="dQ_mean")
    plot(t, t.*fit_dQ, label="Fit")
    legend()
    savefig(string(info.folder, "/Q_profile.pdf"))

end

# Calculate thermal conductivity
function calc_λ(nbins,mbins,dims,dQ_fit)

    T_hot = mbins[1]
    T_cold = mbins[trunc(Int,nbins/2+1)]
    delta_x = dims[1]*0.5
    A = dims[2]*dims[3]

    # λ = dQ /A /((T_hot - T_cold)/delta_x)
    λ = dQ_fit /A /((T_hot-T_cold)/delta_x)

    return λ
end
# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - EvalStructure
# Function to calculate structural properties from position data
# ---
# created by Sebastian Schmitt, 05.04.2020
# ------------------------------------------------------------------------------

## Main function ---------------------------------------------------------------
function eval_structure(dump, info)

    # Calculation of radial distribution function
    calc_rdf(dump,info)

    close("all")
end

## Subfunction -----------------------------------------------------------------
# Function to calculate rdf
function calc_rdf(dump,info)
    # General
    N_bin = info.N_bin
    edges = LinRange(0,info.r_cut,N_bin+1)
    rv = Array((edges[1:end-1] .+ edges[2:end])./2)

    # Get all types and initialize all arrays
    d1 = dump[1]
    types = sort(unique(d1.type))
    combs = []
    for i = 1:length(types), j = i:length(types)
        append!(combs,[[types[i],types[j]]])
    end

    # Calc rdf from all combinations (→ combs)
    N = length(d1.type)
    g_r_all = zeros(N_bin,length(combs)*2)
    Vv = 4/3*π .* (edges[2:end].^3 .- edges[1:end-1].^3)
    N_factor = Float64[]
    for i = 1:length(combs)
        N1 = sum(d1.type .== combs[i][1])
        N2 = sum(d1.type .== combs[i][2])
        if combs[i][1] == combs[i][2]
            push!(N_factor,N1*(N1-1)/2)
        else
            push!(N_factor,N1*N2)
        end
    end
    d1 = nothing

    g_r_all = @distributed (+) for d in dump
        g_r_tmp = zeros(N_bin,length(combs))
        g_r_mol_tmp = zeros(N_bin,length(combs))

        # Box length
        Lx = d.bounds[1,2] - d.bounds[1,1]
        Ly = d.bounds[2,2] - d.bounds[2,1]
        Lz = d.bounds[3,2] - d.bounds[3,1]
        V_box = Lx * Ly * Lz
        ρ_combs = N_factor ./ V_box

        # Loop over all particles
        for i = 1:length(d.type)-1
            # Calculate distances
            rx = abs.(d.x[i] .- d.x[i+1:end]) .% Lx
            rx[rx .> Lx/2] = abs.(rx[rx .> Lx/2] .- Lx)
            ry = abs.(d.y[i] .- d.y[i+1:end]) .% Ly
            ry[ry .> Ly/2] = abs.(ry[ry .> Ly/2] .- Ly)
            rz = abs.(d.z[i] .- d.z[i+1:end]) .% Lz
            rz[rz .> Lz/2] = abs.(rz[rz .> Lz/2] .-Lz)

            r_ij_all = sqrt.( rx.^2 + ry.^2 + rz.^2 )

            # Get all combinations
            combs_i_all = sort([repeat([d.type[i]],N-i) d.type[i+1:end]],dims=2)

            # Exclude all contributions from same molecule
            # → comment out to deactivate {false number of combinations, if activ}
            what_mol = d.molid[i] .!= d.molid[i+1:end]

            combs_i = combs_i_all[what_mol,:]
            r_ij = r_ij_all[what_mol]

            combs_i_mol = combs_i_all[(!).(what_mol),:]
            r_ij_mol = r_ij_all[(!).(what_mol)]

            for q = 1:length(combs)
                # Contributions from atoms from other molecules
                what_q = all(combs_i .== combs[q]',dims=2)[:]
                h = fit(Histogram, r_ij[what_q], edges)
                g_r_tmp[:,q] += h.weights./length(dump)./Vv./ρ_combs[q]

                # Contributions from atoms from the same molecule
                what_q_mol = all(combs_i_mol .== combs[q]',dims=2)[:]
                h_mol = fit(Histogram, r_ij_mol[what_q_mol], edges)
                g_r_mol_tmp[:,q] += h_mol.weights./length(dump)./Vv./ρ_combs[q]
            end
        end
        [g_r_tmp g_r_mol_tmp]
    end
    g_r = g_r_all[:,1:length(combs)]
    g_r_mol = g_r_all[:,length(combs)+1:end]

    # Plot RDFs g(r) for all combinations
    figure()
    for i = 1:length(combs)
        a = combs[i][1]; b = combs[i][2]
        plot(rv,g_r[:,i],label="$a - $b")
    end
    legend()
    if reduced_units
        xlabel(L"r^*")
        ylabel(L"g(r^*)")
    else
        xlabel(L"r / Å")
        ylabel(L"g(r)")
    end
    savefig(string(info.folder,"/fig_rdf.pdf"))

    # Plot RDFs g_mol(r) for all combinations
    figure()
    for i = 1:length(combs)
        a = combs[i][1]; b = combs[i][2]
        plot(rv,g_r_mol[:,i],label="$a - $b")
    end
    legend()
    if reduced_units
        xlabel(L"r^*")
        ylabel(L"g(r^*)")
    else
        xlabel(L"r / Å")
        ylabel(L"g(r)")
    end
    savefig(string(info.folder,"/fig_rdf_mol.pdf"))

    # Output as text file
    filename = string(info.folder,"/rdf.dat")
    header = "r/Å"
    for c in combs
        header = string(header," g_",c[1],c[2])
    end
    for c in combs
        header = string(header," g_mol_",c[1],c[2])
    end
    fID = open(filename,"w")
    println(fID,header)
    writedlm(fID,[rv g_r_all],' ')
    close(fID)
end
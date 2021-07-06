# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - EvalData
# Function to evaluate data from single MD run
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

## Superior Function -----------------------------------------------------------
function EvalSingle(subfolder,inpar)
    # Loading Info File
    moltype, dt, natoms, molmass = load_info(subfolder)

    # Initialization of info structure
    info = info_struct( subfolder,          # info.folder
                        inpar.ensemble,     # info.ensemble
                        inpar.n_equ,        # info.n_equ
                        moltype,            # info.moltype
                        dt,                 # info.dt
                        natoms,             # info.natoms
                        molmass,            # info.molmass
                        inpar.N_bin,        # info.N_bin
                        inpar.r_cut)        # info.r_cut

    # Average Thermodynamic Properties
    T, p, ρ, Etot, Ekin, Epot, c = ave_thermo(info)

    # Evaluate Pressure Data to Calculate Viscosities
    if inpar.mode == "single_run"
        η, η_V = calc_viscosities(info,inpar.corr_length,inpar.span_corr_fun)
    elseif inpar.mode == "tdm"
        η, η_V = calc_viscosities(info)
    end

    # Evaluate Heat Flux Data to Calculate Thermal Conducitvity
    if inpar.mode == "single_run"
        λ = calc_thermalconductivity(info,inpar.corr_length,inpar.span_corr_fun)
    elseif inpar.mode == "tdm"
        λ = calc_thermalconductivity(info)
    end

    # Loading Dump File
    dump = load_dump(info)

    # Evaluate Atoms Positions to calculate Self DIffusivity Coefficient
    D = calc_selfdiffusion(info,dump)

    # Get mole fractions from dump data (first timestep)
    if !(isempty(dump))
        x = get_mole_fraction(info,dump[1])
    else
        x = NaN
    end

    # Structural evaluation
    if inpar.do_structure == 1
        structural = eval_structure(dump,info)
    end

    close("all")

    # Output Results
    OutputResult(results_struct(T, p, ρ, x, Etot, Ekin, Epot, c, η, η_V, D, λ), info.folder)
end

## Function to Average Static Thermodynamic Properties -------------------------
function ave_thermo(info::info_struct)
    # Loading Thermo File
    dat = load_thermo(info)

    what = dat.step .>= info.n_equ

    if info.ensemble == "NVE"
        A = [ones(size(dat.t)) dat.t]
        beta = A\dat.Etot;
        println(info.folder," | ave Etot: ",beta[1]," | slope E_tot: ",beta[2])
    end

    # Temperature
    T = single_dat(mean(dat.T[what]), block_average(dat.T[what])[1], block_average(dat.T[what])[2])
    # Pressure
    if (reduced_units)      factor_p = 1
    elseif !(reduced_units) factor_p = 0.1 end
    p = single_dat(mean(dat.p[what].*factor_p), block_average(dat.p[what])[1].*factor_p, block_average(dat.p[what])[2].*factor_p )
    # Density
    ρ = single_dat(mean(dat.ρ[what]), block_average(dat.ρ[what])[1], block_average(dat.ρ[what])[2])
    # Energies
    Etot = single_dat(mean(dat.Etot[what]), block_average(dat.Etot[what])[1], block_average(dat.Etot[what])[2])
    Ekin = single_dat(mean(dat.Ekin[what]), block_average(dat.Ekin[what])[1], block_average(dat.Ekin[what])[2])
    Epot = single_dat(mean(dat.Epot[what]), block_average(dat.Epot[what])[1], block_average(dat.Epot[what])[2])
    # Heat capacity
    if (reduced_units)      factor_c = 1 / (info.molmass .* info.natoms)
    elseif !(reduced_units) factor_c = eV2J^2 / (kB * (info.molmass*info.natoms/NA/1e3))  end
    c = single_dat((mean(dat.Etot[what]).^2 .- mean(dat.Etot[what].^2.)) ./ T.val^2 * factor_c, NaN, NaN)

    return T, p, ρ, Etot, Ekin, Epot, c
end

## Standard Deviation from Block Average ---------------------------------------
function block_average(x::Array{Float64,1}; M_block=1000)
    # Define maximum Blocklenght and get Number of Blocks
    N_block = trunc(length(x)/M_block)
    # Change M_block if N_block < 20
    if N_block < 20
        N_block = 20
        M_block = trunc(length(x)/N_block)
    end
    # Get Dataset of Blockaverages
    blocks = []
    for k = 1:N_block
        bblock = (round(Int,1+(k-1)*M_block))
        if k == 1
            bblock = 1
        end
        eblock = (round(Int,(k*M_block)))
        if k == N_block
            eblock = round(Int,length(x))
        end
    block = sum(x[bblock:eblock])/M_block
    push!(blocks,block)
    end
    # Get Standard Derivation of Block from Block Averages
    sum_block = 0
    for k = 1:length(blocks)
        sum_block = sum_block + (blocks[k]-mean(x)).^2
    end
    var_block = sum_block/(N_block-1)
    std_block = sqrt(var_block) # Standard Derivation
    err_block = sqrt(var_block/(N_block-1)) # Standard Error
    return std_block, err_block
end


## Functions to apply to atoms positions data
# Function to extract mole fraction from dump data (first timestep)
function get_mole_fraction(info::info_struct, dat::dump_dat)
    Nmol = length(dat.moltype)
    moltype = dat.moltype

    # Calculation of mole fractions on the basis of moltype
    x = Array{single_dat,1}(undef,length(unique(moltype)))
    i = 0
    for imoltype in unique(moltype)
        i += 1
        x[i] = single_dat(sum(moltype .== imoltype) / Nmol,NaN,NaN)
    end
    return x
end

# Function to convert atom data to mol data
function atom2mol(atom)
    N = length(atom)
    mol = Array{Any,1}(undef,N)

    for i = 1:N     # Loop over all timesteps
        molid_atom = atom[i].molid
        mass_atom = atom[i].mass
        if isempty(mass_atom) mass_atom = ones(length(molid_atom)) end
        x_atom = atom[i].x
        y_atom = atom[i].y
        z_atom = atom[i].z

        nmol = maximum(atom[1].molid)
        if !(nmol == length(unique(atom[i].molid)))
            error("Max. molid is not equal to unique vector!")
        end

        molid_mol = Array{Int64,1}(undef,nmol)
        mass_mol = Array{Float64,1}(undef,nmol)
        x_mol = Array{Float64,1}(undef,nmol)
        y_mol = Array{Float64,1}(undef,nmol)
        z_mol = Array{Float64,1}(undef,nmol)
        for imol in unique(molid_atom)  # Loop over all molecules
            # Calculation of the coordinates of COM (center of mass)
            what = molid_atom .== imol
            molid_mol[imol] = imol
            mass_mol[imol] = sum(mass_atom[what])
            x_mol[imol] = (x_atom[what]' * mass_atom[what]) / mass_mol[imol]
            y_mol[imol] = (y_atom[what]' * mass_atom[what]) / mass_mol[imol]
            z_mol[imol] = (z_atom[what]' * mass_atom[what]) / mass_mol[imol]
        end
        mol[i] = dump_dat(atom[i].step,atom[i].t,atom[i].bounds,Int64[],Int64[],molid_mol,atom[1].moltype,mass_mol,x_mol,y_mol,z_mol)
    end

    return mol
end

# Function to extract molecules of given type
function get_mol_by_type(mol,type)
    mol_type = Array{Any,1}(undef,length(mol))
    for i = 1:length(mol)
        imol = mol[i]
        imol_type = dump_dat(imol.step,imol.t,imol.bounds,imol.id,imol.type,imol.molid,imol.moltype,imol.mass,imol.x,imol.y,imol.z)
        what = imol.moltype .== type
        imol_type.molid = imol.molid[what]
        imol_type.moltype = imol.moltype[what]
        imol_type.mass = imol.mass[what]
        imol_type.x = imol.x[what]
        imol_type.y = imol.y[what]
        imol_type.z = imol.z[what]
        mol_type[i] = imol_type
    end
    return mol_type
end
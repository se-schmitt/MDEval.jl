## EvalSingle.jl
# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - EvalData
# Function to evaluate data from single MD run
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------
# [1] Maginn, E. J.; Messerly, R. A.; Carlson, D. J.; Roe, D. R.; Elliott, J. R. Best Practices for Computing Transport Properties 1. Self-Diffusivity and Viscosity from Equilibrium Molecular Dynamics [Article v1.0]. LiveCoMS 2019, 1 (1). https://doi.org/10.33011/livecoms.1.1.6324.

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
                        inpar.n_blocks,     # info.n_blocks
                        inpar.N_bin,        # info.N_bin
                        inpar.r_cut)        # info.r_cut

    # Average Thermodynamic Properties
    T, p, ρ, Etot, Ekin, Epot, c = ave_thermo(info)

    # Calculate box length L_box
    if (reduced_units)
        mass_total = natoms*molmass                         # mass_total = natoms*molmass/NA
        L_box = (mass_total / ρ.val) ^ (1/3)                # [L_box] = 1
    else
        mass_total = natoms*molmass/NA                      # mass_total = natoms*molmass/NA
        L_box = (mass_total / ρ.val * 1e-6) ^ (1/3)         # [L_box] = m
    end

    if inpar.do_transport == 1
        # Evaluate Pressure Data to Calculate Viscosities
        if inpar.mode == "single_run"
            η, η_V = calc_viscosities(info, "single"; mode_acf=inpar.acf_calc_mode, CorrLength=inpar.corr_length, SpanCorrFun=inpar.span_corr_fun, nEvery=inpar.n_every)
        elseif inpar.mode == "tdm"
            η, η_V = calc_viscosities(info, "tdm"; mode_acf=inpar.acf_calc_mode)
        end

        # Evaluate Heat Flux Data to Calculate Thermal Conducitvity
        if inpar.mode == "single_run"
            λ = calc_thermalconductivity(info, "single"; mode_acf=inpar.acf_calc_mode, CorrLength=inpar.corr_length, SpanCorrFun=inpar.span_corr_fun, nEvery=inpar.n_every)
        elseif inpar.mode == "tdm"
            λ = calc_thermalconductivity(info, "tdm"; mode_acf=inpar.acf_calc_mode)
        end
    else
        η = NaN
        η_V = NaN
        λ = NaN
    end

    # Loading Dump File
    dumpexists = false
    if inpar.do_transport == 1 || inpar.do_structure == 1
        dump = load_dump(info)
        dumpexists = true
    end

    if inpar.do_transport == 1
        # Evaluate Atoms Positions to calculate Self DIffusivity Coefficient
        if inpar.mode == "single_run"
            D = calc_selfdiffusion(info, dump; err_mode = "particles")
        elseif inpar.mode == "tdm"
            D = calc_selfdiffusion(info, dump)
        end

        # Finite size correction (for tdm, the correction is applied after averaging the single runs) [1]
        if inpar.mode == "single_run"
            ξ = 2.837298
            if (reduced_units)
                Dcorr = T.val*ξ/(6*π*η.val*L_box)
            else
                Dcorr = kB*T.val*ξ/(6*π*η.val*L_box)
            end
            for i = 1:length(D)
                D[i].val = D[i].val + Dcorr
            end
        end
    else
        D = NaN
    end

    # Get mole fractions from dump data (first timestep)
    if dumpexists
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
    res = results_struct(T, p, ρ, x, Etot, Ekin, Epot, c, η, η_V, D, λ)
    OutputResult(res, info.folder)
end

## Function to Average Static Thermodynamic Properties -------------------------
function ave_thermo(info::info_struct; is_nemd=false)
    # Loading Thermo File
    dat = load_thermo(info; is_nemd=is_nemd)

    what = dat.step .>= info.n_equ

    if info.ensemble == "NVE"
        A = [ones(size(dat.t)) dat.t]
        beta = A\dat.Etot;
        println(info.folder," | ave Etot: ",beta[1]," | slope E_tot: ",beta[2])
    end
    # Temperature
    T_std_err = block_average(dat.T[what],N_blocks=info.n_blocks)
    T = single_dat(mean(dat.T[what]), T_std_err[1], T_std_err[2])
    # Pressure
    if (reduced_units)      factor_p = 1
    elseif !(reduced_units) factor_p = 0.1 end
    p_std_err = block_average(dat.p[what],N_blocks=info.n_blocks)
    p = single_dat(mean(dat.p[what].*factor_p), p_std_err[1].*factor_p, p_std_err[2].*factor_p )
    if is_nemd
        pyz  = single_dat(mean(dat.pyz[what].*factor_p), block_average(dat.pyz[what])[1].*factor_p, block_average(dat.pyz[what])[2].*factor_p)
    else
        pyz = Float64[]
    end
    # Density
    ρ_std_err = block_average(dat.ρ[what],N_blocks=info.n_blocks)
    ρ = single_dat(mean(dat.ρ[what]), ρ_std_err[1], ρ_std_err[2])
    # Energies
    Etot_std_err = block_average(dat.Etot[what],N_blocks=info.n_blocks)
    Etot = single_dat(mean(dat.Etot[what]), Etot_std_err[1], Etot_std_err[2])
    Ekin_std_err = block_average(dat.Ekin[what],N_blocks=info.n_blocks)
    Ekin = single_dat(mean(dat.Ekin[what]), Ekin_std_err[1], Ekin_std_err[2])
    Epot_std_err = block_average(dat.Epot[what],N_blocks=info.n_blocks)
    Epot = single_dat(mean(dat.Epot[what]), Epot_std_err[1], Epot_std_err[2])
    # Heat capacity
    if (reduced_units)      factor_c = 1 / (info.molmass .* info.natoms)
    elseif !(reduced_units) factor_c = eV2J^2 / (kB * (info.molmass*info.natoms/NA/1e3))  end
    c = single_dat((mean(dat.Etot[what]).^2 .- mean(dat.Etot[what].^2.)) ./ T.val^2 * factor_c, NaN, NaN)

    return T, p, ρ, Etot, Ekin, Epot, c, pyz
end

## Standard Deviation from Block Average ---------------------------------------
function block_average(x::Array{Float64,1}; M_block=1000, N_blocks=0)
    if N_blocks == 0
        # Define maximum Blocklenght and get Number of Blocks
        N_blocks = trunc(length(x)/M_block)
        # Change M_block if N_blocks < 20
        if N_blocks < 20
            N_blocks = 20
            M_block = trunc(length(x)/N_blocks)
        end
    elseif N_blocks > 0
        M_block = trunc(length(x)/N_blocks)
    end

    # Get Dataset of Blockaverages
    blocks = []
    for k = 1:N_blocks
        bblock = (round(Int,1+(k-1)*M_block))
        if k == 1
            bblock = 1
        end
        eblock = (round(Int,(k*M_block)))
        if k == N_blocks
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
    var_block = sum_block/(N_blocks-1)
    std_block = sqrt(var_block) # Standard Derivation
    err_block = sqrt(var_block/(N_blocks-1)) # Standard Error
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

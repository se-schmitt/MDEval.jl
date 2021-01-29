# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - EvalData
# Function to evaluate data from single MD run
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

# Main Function
function EvalSingle(info)
    # Loading Info File
    moltype, dt, natoms, molmass = load_info(info.folder)
    info.moltype = moltype
    info.dt = dt
    info.natoms = natoms
    info.molmass = molmass

    # Average Thermodynamic Properties
    T, p, ρ, Etot, Ekin, Epot, c = ave_thermo(info)

    # Evaluate Pressure Data to Calculate Viscosities
    η, η_V = calc_viscosities(info)

    # Evaluate Heat Flux Data to Calculate Thermal Conducitvity
    λ = calc_thermalconductivity(info)

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
    close("all")

    # Output Results
    OutputResult(results_struct(T, p, ρ, x, Etot, Ekin, Epot, c, η, η_V, D, λ), info.folder)
end

# Subfunctions
# Function to Average Thermodynamic Properties
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

# Evaluate Pressure Data to calculate viscosities
function calc_viscosities(info::info_struct)
    # Loading Pressure File
    dat = load_pressure(info)

    if (dat isa pressure_dat)
        what = (dat.step .>= info.n_equ)
        t = dat.t[what]

        # Calculate autocorrelation
        P_mat = [dat.pxy[what],
                 dat.pxz[what],
                 dat.pyz[what],
                 dat.pxx[what] .- dat.pyy[what],
                 dat.pyy[what] .- dat.pzz[what],
                 dat.pxx[what] .- dat.pzz[what]]
        acf_P = pmap(calc_acf,P_mat)

        # Unit conversion and GK factor
        metal2Pas = 1e-32
        factor = metal2Pas .* dat.V[what] ./ (kB .* dat.T[what])
        if (reduced_units) factor = dat.V[what] ./ (dat.T[what]) end

        # Integration of acf
        acf_η = (acf_P[1] .+ acf_P[2] .+ acf_P[3])./6 .+
                (acf_P[4] .+ acf_P[5] .+ acf_P[6])./24
        η_t = cumtrapz(t,acf_η).*factor

        # Bulk viscosity
        P_diag = [dat.pxx[what], dat.pyy[what], dat.pzz[what]]
        acf_η_V = mean(pmap(calc_acf,P_diag))
        η_V_t = cumtrapz(t,acf_η_V).*factor

        # Write data to file
        line1 = string("# Created by MD - Bulk Evaluation, Folder: ", info.folder)
        line2 = "# t[ps] η(t)[Pa*s] ACF_η[Pa^2] η_V(t)[Pa*s] ACF_η_V[Pa^2]"
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

# Evaluate Atoms Positions to calculate Self DIffusivity Coefficient
function calc_selfdiffusion(info::info_struct, dat::Array{Any,1})
    if (!isempty(dat))
        N_moltype = maximum(dat[1].moltype)
        D = Array{single_dat,1}(undef,N_moltype)

        # Convert atom to molecule coordinates
        mol = atom2mol(dat)

        # Initialize figure and msd_t_all variable
        colors = ["b","g","r","c","m","y"]
        figure()
        msd_t_all = Array{Float64,2}(undef,0,N_moltype)

        for i = 1:N_moltype
            # Get all molecules of type 'moltype'
            mol_i = get_mol_by_type(mol,i)

            # Calculation of Mean Square Displacement
            msd_t, t = calc_msd(mol_i,info)

            # Segmenting the msd and calculation of log-log slope
            kmax_eval = 0.6
            δ_tol = 0.075
            δ = 0
            N_eval = ceil(Int64,kmax_eval*length(t))
            L_eval = 50
            N = N_eval
            while δ < δ_tol && N > L_eval
                N = N-L_eval
                what = Array(N+1:N+L_eval)
                beta = hcat(ones(L_eval),log10.(t[what]))\log10.(msd_t[what])
                δ = abs(beta[2]-1)
            end

            what_eval = Int64.(N+1:N_eval)
            if length(what_eval) < 150
                what_eval = Int64.(ceil(0.25*length(t)):ceil(0.75*length(t)))
            end

            # Calculation of D by determination of slope of msd
            beta = hcat(ones(length(what_eval)),t[what_eval]) \ msd_t[what_eval]
            factor_unit = 1e-8
            if (reduced_units) factor_unit = 1 end
            D[i] = single_dat(beta[2]/6*factor_unit,NaN,NaN)

            # Plot MSD of molecule i
            t1 = round(Int64,t[what_eval[1]])
            t2 = round(Int64,t[what_eval[end]])

            plot(t, msd_t, linewidth=1, color=colors[i], label=string("Mol. ",i," (",t1," ps - ",t2," ps (",length(what_eval)," steps))"))
            plot(t,beta[1] .+ beta[2].*t, linestyle="-", linewidth=1, color=colors[i])
            plot([t1,t1],[0.1,msd_t[what_eval[1]]], color="black", linestyle=":", linewidth=1)
            plot([t2,t2],[0.1,msd_t[what_eval[end]]], color="black", linestyle=":", linewidth=1)

            # Save MSD in array to save dlm file
            if i == 1
                msd_t_all = t
            end
            msd_t_all = hcat(msd_t_all, msd_t)
        end

        # Save msd as data file
        line1 = string("# Created by MD - Bulk Evaluation, Folder: ", info.folder)
        line2 = string("# t[ps]"," msd[Å]"^N_moltype)
        header = string(line1,"\n",line2)
        file = string(info.folder,"/msd.dat")
        fID = open(file,"w"); println(fID,header)
        writedlm(fID, msd_t_all, " ")
        close(fID)

        # Formatting figure
        if !(reduced_units) xlabel(L"t / ps") elseif reduced_units xlabel(L"t*") end
        if !(reduced_units) ylabel(L"MSD / Å$^2$") elseif reduced_units ylabel(L"MSD*") end
        title("Mean square displacement (MSD)")
        legend()
        tight_layout()
        savefig(string(info.folder,"/fig_msd(t).pdf"))
        close()
    else
        D = single_dat(NaN,NaN,NaN)
    end

    return D
end

# Evaluate Heat Flux Data to Calculate Thermal Conducitvity
function calc_thermalconductivity(info::info_struct)
    # Loading Heat Flux File
    dat = load_heatflux(info)

    if (dat isa heat_dat)
        what = (dat.step .>= info.n_equ)
        t = dat.t[what]

        # Calculate autocorrelation
        J_mat = [dat.jx[what]./dat.V[what],
                 dat.jy[what]./dat.V[what],
                 dat.jz[what]./dat.V[what]]
        acf_J = pmap(calc_acf,J_mat)

        # Unit conversion and GK factor
        metal2WmK = 2.566969967e-16
        factor = metal2WmK .* dat.V[what] ./ (3 .* kB .* dat.T[what].^2)
        if (reduced_units) factor = dat.V[what] ./ (3 .* dat.T[what].^2) end

        # Integration of acf
        acf_λ = (acf_J[1] .+ acf_J[2] .+ acf_J[3])
        λ_t = cumtrapz(t,acf_λ).*factor

        # Write data to file
        line1 = string("# Created by MD - Bulk Evaluation, Folder: ", info.folder)
        line2 = "# t[ps] λ(t)[W/(m*K)] ACF_λ[eV^2/(Å^4*ps)]"
        header = string(line1,"\n",line2)
        if (reduced_units) line2 = "# t* λ*(t) ACF_λ*" end
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
        savefig(string(info.folder,"/fig_lambda_acf(t).pdf"))
        close()
        figure()
        plot(t[1:teve:end],λ_t[1:teve:end], linewidth=0.1)
        if !(reduced_units) xlabel(L"t / ps") elseif reduced_units xlabel(L"t*") end
        if !(reduced_units) ylabel(L"λ / (W/(m*K)") elseif reduced_units ylabel(L"λ*") end
        savefig(string(info.folder,"/fig_lambda(t).pdf"))
        close()

        λ = single_dat(λ_t[end],NaN,NaN)
    else
        λ = single_dat(NaN,NaN,NaN)
    end

    return λ
end

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

# Help functions
# Function to integrate by trapezodial rule
function cumtrapz(x::Array{Float64},y::Array{Float64})
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

# Function to calculate autocorrelation function (usable parallel)
@everywhere function calc_acf(in::Array{Float64,1})
    L = length(in)
    in0 = vcat(in,zeros(L))
    out = real(ifft(fft(in0) .* conj(fft(in0))))[1:L]
    out = out ./ (L:-1:1)
    return out
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

# Function to calculate mean square displacement
function calc_msd(dat::Array{Any,1},info)
    n = length(dat[1].molid)    # No. molecules
    N = length(dat)             # No. timesteps

    t = Float64[]
    for tmp in dat append!(t,tmp.t) end
    x = fill(Float64[],n)
    y = fill(Float64[],n)
    z = fill(Float64[],n)
    i = 0
    for imol in dat[1].molid
        i += 1
        xmol = Float64[]
        ymol = Float64[]
        zmol = Float64[]
        for tmp in dat
            what = tmp.molid .== imol
            append!(xmol,tmp.x[what])
            append!(ymol,tmp.y[what])
            append!(zmol,tmp.z[what])
        end
        x[i] = xmol; y[i] = ymol; z[i] = zmol
    end

    msd_t = zeros(N)

    for i = 1:n
        r2 = x[i] .^ 2 + y[i] .^ 2 + z[i] .^ 2

        s_AB = sum(pmap(calc_acf,[x[i],y[i],z[i]]))

        msd_i = Array{Float64,1}(undef,N)
        sumsq = []
        for im = 1:N
            m = im-1
            if (im == 1) sumsq = 2*sum(r2)
            else         sumsq = sumsq - r2[im-1] - r2[N-m] end
            msd_i[im] = sumsq / (N - m) - 2*s_AB[im]
        end

        msd_t += msd_i ./n
    end
    msd_t = abs.(msd_t)

    return msd_t, t
end

# Standard Deviation from Block Average
function block_average(x::Array{Float64,1})
    # Define maximum Blocklenght and get Number of Blocks
    M_block = 1000
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

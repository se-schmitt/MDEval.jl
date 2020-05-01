# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - EvalData
# Function to evaluate data from single MD run
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

# Main Function
function EvalData(info)
    # Loading Info File
    moltype, dt, natoms, molmass = load_info(info.folder)
    info.moltype = moltype
    info.dt = dt
    info.natoms = natoms
    info.molmass = molmass

    # Average Thermodynamic Properties
    T, p, ρ, Etot, Ekin, Epot = ave_thermo(info)

    # Evaluate Pressure Data to Calculate Viscosities
    η, η_V = calc_viscosities(info)

    # Evaluate Atoms Positions to calculate Self DIffusivity Coefficient
    D = calc_selfdiffusion(info)

    # Evaluate Heat Flux Data to Calculate Thermal Conducitvity
    λ = calc_thermalconductivity(info)

    # Output Results
    OutputResult(results_struct(T, p, ρ, Etot, Ekin, Epot, η, η_V, D, λ), info.folder)
end

# Subfunctions
# Function to Average Thermodynamic Properties
function ave_thermo(info::info_struct)
    # Loading Thermo File
    dat = load_thermo(info)

    what = dat.step .>= info.n_equ#

    if info.ensemble == "NVE"
        A = [ones(size(dat.t)) dat.t]
        beta = A\dat.Etot;
        println(info.folder," | ave Etot: ",beta[1]," | slope E_tot: ",beta[2])
    end

    T = single_dat(mean(dat.T[what]), std(dat.T[what]), NaN)
    p = single_dat(mean(dat.p[what].*1e-1), std(dat.p[what]).*1e-1, NaN)
    ρ = single_dat(mean(dat.ρ[what]), std(dat.ρ[what]), NaN)
    Etot = single_dat(mean(dat.Etot[what]), std(dat.Etot[what]), NaN)
    Ekin = single_dat(mean(dat.Ekin[what]), std(dat.Ekin[what]), NaN)
    Epot = single_dat(mean(dat.Epot[what]), std(dat.Epot[what]), NaN)

    return T, p, ρ, Etot, Ekin, Epot
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
        header = string(line1,"\n",line2)
        file = string(info.folder,"/viscosity.dat")
        fID = open(file,"w"); println(fID,header)
        writedlm(fID, hcat(dat.t[what],η_t,acf_η,η_V_t,acf_η_V), " ")
        close(fID)

        # Plots
        teve = 100        # Just take data every 'teve' timesteps
        plt = plot(t[2:teve:end],abs.(acf_η[2:teve:end]), dpi=400, legend=false, xscale=:log10, yscale=:log10)
        xlabel!("t / ps"); ylabel!("|ACF| / Pa^2")
        png(plt,string(info.folder,"/fig_eta-acf(t).png"))
        plt = plot(t[1:teve:end],η_t[1:teve:end], dpi=400, legend=false)
        xlabel!("t / ps"); ylabel!("η / (Pa*s)")
        png(plt,string(info.folder,"/fig_eta(t).png"))

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
function calc_selfdiffusion(info::info_struct)
    # Loading Dump File
    dat = load_dump(info)

    if (!isempty(dat))
        # Convert atom to molecule coordinates
        mol = atom2mol(dat)

        # Calculation of Mean Square Displacement
        msd_t, t = calc_msd(mol,info)

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
        D = single_dat(beta[2]/6*factor_unit,NaN,NaN)

        t1 = t[what_eval[1]]
        t2 = t[what_eval[end]]

        # Plots of msd
        infostr = string("MSD 'Evaluation window':\n",t1," ps - ",t2," ps (",length(what_eval)," steps))")
        plt1 = plot(t,msd_t, dpi=400, legend=false, title=infostr, titlefont=font(10))
        plot!(t,beta[1] .+ beta[2].*t,linestyle=:dash)
        plot!([t1,t1],[0.1,msd_t[what_eval[1]]],linestyle=:dot)
        plot!([t2,t2],[0.1,msd_t[what_eval[end]]],linestyle=:dot)
        png(plt1,string(info.folder,"/fig_msd(t).png"))
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

        # Integration of acf
        acf_λ = (acf_J[1] .+ acf_J[2] .+ acf_J[3])
        λ_t = cumtrapz(t,acf_λ).*factor

        # Write data to file
        line1 = string("# Created by MD - Bulk Evaluation, Folder: ", info.folder)
        line2 = "# t[ps] λ(t)[W/(m*K)] ACF_λ[eV^2/(Å^4*ps)]"
        header = string(line1,"\n",line2)
        file = string(info.folder,"/thermalconductivity.dat")
        fID = open(file,"w"); println(fID,header)
        writedlm(fID, hcat(dat.t[what],λ_t,acf_λ), " ")
        close(fID)

        # Plots
        teve = 100        # Just take data every 'teve' timesteps
        plt = plot(t[2:teve:end],abs.(acf_λ[2:teve:end]), dpi=400, legend=false, xscale=:log10, yscale=:log10)
        xlabel!("t / ps"); ylabel!("|ACF| / eV^2/(Å^4*ps)")
        png(plt,string(info.folder,"/fig_lambda_acf(t).png"))
        plt = plot(t[1:teve:end],λ_t[1:teve:end], dpi=400, legend=false)
        xlabel!("t / ps"); ylabel!("λ / (W/(m*K))")
        png(plt,string(info.folder,"/fig_lambda(t).png"))

        λ = single_dat(λ_t[end],NaN,NaN)
    else
        λ = single_dat(NaN,NaN,NaN)
    end

    return λ
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
        mass_mol = Array{Int64,1}(undef,nmol)
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
        mol[i] = dump_dat(atom[i].step,atom[i].t,atom[i].bounds,[],molid_mol,mass_mol,x_mol,y_mol,z_mol)
    end

    return mol
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
    for imol = 1:n
        xmol = Float64[]
        ymol = Float64[]
        zmol = Float64[]
        for tmp in dat
            what = tmp.molid .== imol
            append!(xmol,tmp.x[what])
            append!(ymol,tmp.y[what])
            append!(zmol,tmp.z[what])
        end
        x[imol] = xmol; y[imol] = ymol; z[imol] = zmol
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

    return msd_t, t
end

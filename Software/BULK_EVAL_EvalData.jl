# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - EvalData
# Function to evaluate data from single MD run
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

# Main Function
function EvalData(DATA, info)

    thermodat, pdat, posdat, jdat = DATA2dats(DATA)

    # Average Thermodynamic Properties
    T, p, ρ, Etot, Ekin, Epot = ave_thermo(thermodat, info)


    # Evaluate Pressure Data to Calculate Viscosities
    if (pdat isa pressure_dat) η, η_V = calc_viscosities(pdat, info)
    else η = single_dat([],[],[]); η_V = single_dat([],[],[]) end

    # Evaluate Atoms Positions to calculate Self DIffusivity Coefficient
    if (!isempty(posdat)) D = calc_selfdiffusion(posdat, info)
    else D = single_dat([],[],[]) end

    # Evaluate Heat Flux Data to Calculate Thermal Conducitvity
    if (jdat isa heat_dat) λ = calc_thermalconductivity(jdat, info)
    else λ = single_dat([],[],[]) end

    RES = results_struct(T, p, ρ, Etot, Ekin, Epot, η, η_V, D, λ)
    return RES
end

# Subfunctions
# Function to Average Thermodynamic Properties
function ave_thermo(dat, info)
    what = dat.step .>= info.n_equ

    T = single_dat(mean(dat.T[what]), std(dat.T[what]), [])
    p = single_dat(mean(dat.p[what].*1e-1), std(dat.p[what]).*1e-1, [])
    ρ = single_dat(mean(dat.ρ[what]), std(dat.ρ[what]), [])
    Etot = single_dat(mean(dat.Etot[what]), std(dat.Etot[what]), [])
    Ekin = single_dat(mean(dat.Ekin[what]), std(dat.Ekin[what]), [])
    Epot = single_dat(mean(dat.Epot[what]), std(dat.Epot[what]), [])

    return T, p, ρ, Etot, Ekin, Epot
end

# Evaluate Pressure Data to calculate viscosities
function calc_viscosities(dat, info)            # dat => pdat
    what = (dat.step .>= info.n_equ) .& (dat.step .<= 2e6)
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
    png(plt,string(info.folder,"/fig_acf(t).png"))
    plt = plot(t[1:teve:end],η_t[1:teve:end], dpi=400, legend=false)
    xlabel!("t / ps"); ylabel!("η / (Pa*s)")
    png(plt,string(info.folder,"/fig_eta(t).png"))

    # Save Results
    η = single_dat(η_t[end], [], [])
    η_V = single_dat(η_V_t[end], [], [])
    return η, η_V
end

# Evaluate Atoms Positions to calculate Self DIffusivity Coefficient
function calc_selfdiffusion(dat, info)          # dat => posdat
    return single_dat([],[],[])
end

# Evaluate Heat Flux Data to Calculate Thermal Conducitvity
function calc_thermalconductivity(dat, info)    # dat => jdat
    return single_dat([],[],[])
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

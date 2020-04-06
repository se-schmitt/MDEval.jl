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
    else η = single_dat([],[]); η_V = single_dat([],[]) end

    # Evaluate Atoms Positions to calculate Self DIffusivity Coefficient
    if (!isempty(posdat)) D = calc_selfdiffusion(posdat, info)
    else D = single_dat([],[]) end

    # Evaluate Heat Flux Data to Calculate Thermal Conducitvity
    if (jdat isa heat_dat) λ = calc_thermalconductivity(jdat, info)
    else λ = single_dat([],[]) end

    RES = results_struct(T, p, ρ, Etot, Ekin, Epot, η, η_V, D, λ)
    return RES
end

# Subfunctions
# Function to Average Thermodynamic Properties
function ave_thermo(dat, info)
    what = dat.step .>= info.n_equ

    T = single_dat(mean(dat.T[what]), std(dat.T[what]))
    p = single_dat(mean(dat.p[what].*1e-1), std(dat.p[what]).*1e-1)
    ρ = single_dat(mean(dat.ρ[what]), std(dat.ρ[what]))
    Etot = single_dat(mean(dat.Etot[what]), std(dat.Etot[what]))
    Ekin = single_dat(mean(dat.Ekin[what]), std(dat.Ekin[what]))
    Epot = single_dat(mean(dat.Epot[what]), std(dat.Epot[what]))

    return T, p, ρ, Etot, Ekin, Epot
end

# Evaluate Pressure Data to calculate viscosities
function calc_viscosities(dat, info)            # dat => pdat
    what = dat.step .>= info.n_equ

    # Create pressure tensor
    tr_p = (dat.pxx .+ dat.pyy .+ dat.pzz)./3
    P = hcat(dat.pxx,dat.pyy,dat.pzz,dat.pxy,dat.pxz,dat.pyz)[what,:]

    # Calculate autocorrelation
    lags = [0:size(P,1)-1;]
    acf_P = autocor(P, lags; demean=false)

    # Integration of acf
    metal2Pas = 1e-32
    factor = metal2Pas .* dat.V[what] ./ (kB .* dat.T[what])
    η_tensor = trapz(dat.t[what], acf_P, 1) .* repeat(factor,1,6)

    # Shear viscosity
    # ave_method: "off-diag" - just off-diagonal components used
    #             "all"      - all components used (ff. Daivies and Evans 1994)
    ave_method = "off-diag"

    if ave_method == "off-diag"
        η_t = mean(η_tensor[:,4:6],dims=2)
        acf_η = mean(acf_P[:,4:6],dims=2).*1e10
    end

    # Bulk viscosity
    η_V_t = mean(η_tensor[:,1:3],dims=2)
    acf_η_V = mean(acf_P[:,1:3],dims=2).*1e10

    # Write data to file
    header = "# t[ps] ACF_η[Pa^2] η(t)[Pa*s] ACF_η_V[Pa^2] η_V(t)[Pa*s]"
    file = string(info.folder,"/viscosity.dat")
    fID = open(file,"w"); println(fID,header)
    writedlm(fID, hcat(dat.t[what],acf_η,η_t,acf_η_V,η_V_t), " ")
    close(fID)

    # Plots
    plt = plot(dat.t[what],abs.(acf_η), dpi=400, legend=false, yscale=:log)
    xlabel!("t / ps"); ylabel!("|ACF| / Pa^2")
    png(plt,string(info.folder,"/fig_acf(t).png"))
    plt = plot(dat.t[what],η_t, dpi=400, legend=false)
    xlabel!("t / ps"); ylabel!("η / (Pa*s)")
    png(plt,string(info.folder,"/fig_eta(t).png"))

    # Save Results
    η = single_dat(η_t[end],[])
    η_V = single_dat(η_V_t[end],[])
    return η, η_V
end

# Evaluate Atoms Positions to calculate Self DIffusivity Coefficient
function calc_selfdiffusion(dat, info)          # dat => posdat
    return single_dat([],[])
end

# Evaluate Heat Flux Data to Calculate Thermal Conducitvity
function calc_thermalconductivity(dat, info)    # dat => jdat
    return single_dat([],[])
end

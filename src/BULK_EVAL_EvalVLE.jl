# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - EvalVLE
# Function to evaluate VLE simulations
# ---
# created by Sebastian Schmitt, 03.09.2020
# ------------------------------------------------------------------------------

function EvalVLE(info)
    # Get subfolder
    subfolder = get_subfolder(info.folder)

    # Evluation of subfolders
    for f in subfolder
        # Loading Info File
        moltype, dt, natoms, molmass = load_info(f)
        info.moltype = moltype
        info.dt = dt
        info.natoms = natoms
        info.molmass = molmass

        # Evaluation of thermo files
        eval_thermo_vle(f,info)

        # Evaluation of profile files
        eval_profiles(f,info)
        println(string("Subfolder ",findfirst(f .== subfolder)," / ",length(subfolder)," DONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS")))
    end

    # Plot of the vapour and liquid densities
    T = Float64[]
    ρ = Float64[]
    ρv = Float64[]
    ρl = Float64[]
    px = Float64[]
    py = Float64[]
    pz = Float64[]
    D = Float64[]

    for f in subfolder
        file_res = string(f,"/result.dat")
        fID = open(file_res,"r")
        txt = readlines(fID)
        close(fID)
        for line in txt
            if line[1] != "#"
                pos0 = findfirst(isequal(':'),line)
                pos1 = findfirst(isequal(' '),line)+1
                pos2 = findlast(isequal(','),line)-1

                if line[1:pos0] == "T:"         append!(T,parse(Float64,strip(line[pos1:pos2])))
                elseif line[1:pos0] == "ρ:"     append!(ρ,parse(Float64,strip(line[pos1:pos2])))
                elseif line[1:pos0] == "ρv:"    append!(ρv,parse(Float64,strip(line[pos1:pos2])))
                elseif line[1:pos0] == "ρl:"    append!(ρl,parse(Float64,strip(line[pos1:pos2])))
                elseif line[1:pos0] == "px:"    append!(px,parse(Float64,strip(line[pos1:pos2])))
                elseif line[1:pos0] == "py:"    append!(py,parse(Float64,strip(line[pos1:pos2])))
                elseif line[1:pos0] == "pz:"    append!(pz,parse(Float64,strip(line[pos1:pos2])))
                elseif line[1:pos0] == "D:"     append!(D,parse(Float64,strip(line[pos1:pos2]))) end
            end
        end
    end

    header = "T/K ρ/(g/ml) ρ_v/(g/ml) ρ_l/(g/ml) px/MPa py/MPa pz/MPa D/Å"
    ind = sortperm(T)
    mat = hcat(T,ρ,ρv,ρl,px,py,pz,D)[ind,:]

    # Data file
    fID = open(string(info.folder,"/results.dat"),"w")
    println(fID,header)
    writedlm(fID,mat," ")
    close(fID)

    # Plot
    figure()
    for i = 1:length(T)
        plot([ρv[i],ρl[i]],[T[i],T[i]],"sb:",linewidth=0.5)
    end
    xlabel(L"\rho_{m} / (g/ml)")
    ylabel(L"T / K")
    title("Vapour-liquid equilibrium")
    savefig(string(info.folder,"/Fig_rho-T_VLE.pdf"))

    close("all")
end

## Subfunctions
# Function to evaluate thermo file
function eval_thermo_vle(folder,info)
    # Load data
    dat = load_thermo_vle(folder, info)
    dat.px = dat.px .* 0.1
    dat.py = dat.py .* 0.1
    dat.pz = dat.pz .* 0.1

    # Output data
    OutputResult_VLE(dat, folder)
end

# Function to evaluate profiles of density, temperature, stress
function eval_profiles(folder,info)
    # Load profile file
    dat = []
    ts_add = 0
    for i = 1:10
        filename = string(folder,"/vle.2phase.",i,".profile")
        if isfile(filename)
            dat = read_profile1D(filename,dat,ts_add)
            ts_add = dat.timestep[end]
        end
    end

    # Calculation of means
    x = mean(dat.x,dims=1)[:]
    x = x .- minimum(x)
    Ncount = mean(dat.Ncount,dims=1)[:]
    ρn = mean(dat.ρn,dims=1)[:]
    ρm = mean(dat.ρm,dims=1)[:]
    T = sum(dat.T .* dat.Ncount,dims=1)[:] ./ sum(dat.Ncount,dims=1)[:]
    pxx = sum(dat.pxx .* dat.Ncount,dims=1)[:] ./ sum(dat.Ncount,dims=1)[:] .*0.1
    pyy = sum(dat.pyy .* dat.Ncount,dims=1)[:] ./ sum(dat.Ncount,dims=1)[:] .*0.1
    pzz = sum(dat.pzz .* dat.Ncount,dims=1)[:] ./ sum(dat.Ncount,dims=1)[:] .*0.1
    pxy = sum(dat.pxy .* dat.Ncount,dims=1)[:] ./ sum(dat.Ncount,dims=1)[:] .*0.1
    pxz = sum(dat.pxz .* dat.Ncount,dims=1)[:] ./ sum(dat.Ncount,dims=1)[:] .*0.1
    pyz = sum(dat.pyz .* dat.Ncount,dims=1)[:] ./ sum(dat.Ncount,dims=1)[:] .*0.1

    # pxx = mean(dat.pxx,dims=1)[:].*0.1
    # pyy = mean(dat.pyy,dims=1)[:].*0.1
    # pzz = mean(dat.pzz,dims=1)[:].*0.1
    # pxy = mean(dat.pxy,dims=1)[:].*0.1
    # pxz = mean(dat.pxz,dims=1)[:].*0.1
    # pyz = mean(dat.pyz,dims=1)[:].*0.1

    # Get position of interfaces
    ρm_ = (maximum(ρm) - minimum(ρm))/2 + minimum(ρm)
    x_flu = x[ρm .> ρm_]
    ind = sortperm(ρm)
    ρm_sorted = ρm[ind]
    x_sorted = x[ind]
    what = x_sorted.<mean(x_flu)
    fun_int1 = Spline1D(ρm_sorted[what],x_sorted[what],k=1)
    x_interface1 = fun_int1(ρm_)
    what = x_sorted.>mean(x_flu)
    fun_int2 = Spline1D(ρm_sorted[what],x_sorted[what],k=1)
    x_interface2 = fun_int2(ρm_)

    x_flu_mid = (x_interface2+x_interface1)/2
    part1 = x .< x_flu_mid
    part2 = x .>= x_flu_mid

    dx = zeros(size(x))
    dx[part1] = x[part1] .- x_interface1
    dx[part2] = -1 .* (x[part2] .- x_interface2)

    # Fitting of density profile
    # Parameter: p[1] = ρ_l, p[2] = ρ_v, p[3] = D
    # Classical
    @. fun_profile(x, p) = 0.5 .*((p[1].+p[2]) .- (p[1].-p[2]) .* tanh.(2 .* (x.-p[4]) ./ p[3]))

    p0 = [maximum(ρm),minimum(ρm),10,0]
    fit = curve_fit(fun_profile, dx, ρm, p0)

    if fit.param[1] > fit.param[2]
        ρ_l = fit.param[1]
        ρ_v = fit.param[2]
        case = 1
    else
        ρ_v = fit.param[1]
        ρ_l = fit.param[2]
        case = 2
    end
    D = fit.param[3]

    if case == 1
        what_l = dx .< -2*D
        what_v = dx .> 2*D
    elseif case == 2
        what_l = dx .> 2*D
        what_v = dx .< -2*D
    end

    p_l = (pxx[what_l] .+ pyy[what_l] .+ pzz[what_l])' * Ncount[what_l] ./(sum(Ncount[what_l]).*3)
    p_v = (pxx[what_v] .+ pyy[what_v] .+ pzz[what_v])' * Ncount[what_v] ./(sum(Ncount[what_v]).*3)
    p_x = pxx'*Ncount/sum(Ncount)
    Tm = T'*Ncount/sum(Ncount)

    # Data file
    fID = open(string(folder,"/result.dat"),"a")
    print_prop(fID, single_dat(Tm,NaN,NaN), "Tpro")
    print_prop(fID, single_dat(ρ_l,NaN,NaN), "ρl")
    print_prop(fID, single_dat(ρ_v,NaN,NaN), "ρv")
    print_prop(fID, single_dat(D,NaN,NaN), "D")
    close(fID)

    # Figure
    figure()
    plot(dx,ρm,"*")
    plot(dx[what_l],ρm[what_l],"o",fillstyle="none")
    plot(dx[what_v],ρm[what_v],"s",fillstyle="none")
    xfit = Array(minimum(dx):0.001:maximum(dx))
    plot(xfit,fun_profile(xfit, fit.param))
    legend([L"\rho_{m}: \mathrm{all}",L"\rho_{m}: \mathrm{liquid-bulk}",L"\rho_{m}: \mathrm{vapor-bulk}",L"\rho_{m}: \mathrm{fit}"],loc="lower left")
    title("Mass density + Fit")
    xlabel(L"\Delta{}x / \mathrm{\AA}")
    ylabel(L"\rho_\mathrm{m} / \mathrm{(g/ml)}")
    grid();    tight_layout()
    savefig(string(folder,"/Fig_mdensity_fit.pdf"))

    figure()
    plot(dx,T)
    title("Temperature")
    xlabel(L"\Delta{}x / \mathrm{\AA}")
    ylabel(L"T / \mathrm{K}")
    grid();    tight_layout()
    savefig(string(folder,"/Fig_temperature.pdf"))

    figure()
    plot(x,ρn,".-")
    title("Number density")
    xlabel(L"x / \mathrm{\AA}")
    ylabel(L"\rho_\mathrm{n} / \mathrm{\AA^{3}}")
    grid();    tight_layout()
    savefig(string(folder,"/Fig_ndensity.pdf"))

    figure()
    plot(x,ρm,".-")
    plot([x_interface1,x_interface2],[ρm_,ρm_],"*-")
    title("Mass density")
    xlabel(L"x / \mathrm{\AA}")
    ylabel(L"\rho_\mathrm{m} / \mathrm{(g/ml)}")
    grid();    tight_layout()
    savefig(string(folder,"/Fig_mdensity.pdf"))

    figure()
    plot(dx,pxx);    plot(dx,pyy);    plot(dx,pzz)
    plot(dx,pxy);    plot(dx,pxz);    plot(dx,pyz)
    title("Stress")
    xlabel(L"\Delta{}x / \mathrm{\AA}")
    ylabel(L"p / \mathrm{MPa}")
    legend([L"p_{xx}",L"p_{yy}",L"p_{zz}",L"p_{xy}",L"p_{xz}",L"p_{yz}"],loc="lower left")
    grid();    tight_layout()
    savefig(string(folder,"/Fig_stress.pdf"))

    figure()
    contourf(x, dat.timestep, dat.ρm)
    xlabel(L"x / \mathrm{\AA}")
    ylabel(L"N_\mathrm{step}")
    f = matplotlib.ticker.FormatStrFormatter("%1.1e")
    ax = gca()
    ax.yaxis.set_major_formatter(f)
    tight_layout()
    savefig(string(folder,"/Fig_density_f(ts).pdf"))

    close("all")
end

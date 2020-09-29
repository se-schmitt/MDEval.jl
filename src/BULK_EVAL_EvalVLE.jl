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
        eval_profiles(f)
        println(string("Subfolder ",findfirst(f .== subfolder)," / ",length(subfolder)," DONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS")))
    end

    # Plot of the vapour and liquid densities
    T = Float64[]
    ρ_v = Float64[]
    ρ_l = Float64[]
    p_v = Float64[]
    p_l = Float64[]
    D = Float64[]

    for f in subfolder
        file_res = string(f,"/result.dat")
        fID = open(file_res,"r")
        txt = readlines(fID)
        close(fID)
        for line in txt
            pos0 = findfirst(isequal(':'),line)-1
            pos1 = findfirst(isequal(' '),line)+1
            pos2 = findlast(isequal(' '),line)-1
            if line[1:pos0] == "T"          append!(T,parse(Float64,line[pos1:pos2]))
            elseif line[1:pos0] == "ρ_v"    append!(ρ_v,parse(Float64,line[pos1:pos2]))
            elseif line[1:pos0] == "ρ_l"    append!(ρ_l,parse(Float64,line[pos1:pos2]))
            elseif line[1:pos0] == "p_v"    append!(p_v,parse(Float64,line[pos1:pos2]))
            elseif line[1:pos0] == "p_l"    append!(p_l,parse(Float64,line[pos1:pos2]))
            elseif line[1:pos0] == "D"      append!(D,parse(Float64,line[pos1:pos2])) end
        end
    end

    header = "T/K ρ_v/(g/ml) ρ_l/(g/ml) p_v/MPa p_l/MPa D/Å"
    ind = sortperm(T)
    mat = hcat(T,ρ_v,ρ_l,p_v,p_l,D)[ind,:]

    # Data file
    fID = open(string(info.folder,"/results.dat"),"w")
    println(fID,header)
    writedlm(fID,mat," ")
    close(fID)

    # Plot
    figure()
    for i = 1:length(T)
        plot([ρ_v[i],ρ_l[i]],[T[i],T[i]],"sb:",linewidth=0.5)
    end
    xlabel(L"\rho_{m} / (g/ml)")
    ylabel(L"T / K")
    title("Vapour-liquid equilibrium")
    savefig(string(info.folder,"/Fig_rho-T_VLE.pdf"))

    close("all")
end

## Subfunctions
# Function to evaluate profiles of density, temperature, stress
function eval_profiles(folder)
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
    ρn = mean(dat.ρn,dims=1)[:]
    ρm = mean(dat.ρm,dims=1)[:]
    T = mean(dat.T,dims=1)[:]
    pxx = mean(dat.pxx,dims=1)[:]
    pyy = mean(dat.pyy,dims=1)[:]
    pzz = mean(dat.pzz,dims=1)[:]
    pxy = mean(dat.pxy,dims=1)[:]
    pxz = mean(dat.pxz,dims=1)[:]
    pyz = mean(dat.pyz,dims=1)[:]

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

    p_l = mean(pxx[what_l] .+ pyy[what_l] .+ pzz[what_l])/3
    p_v = mean(pxx[what_v] .+ pyy[what_v] .+ pzz[what_v])/3

    # Data file
    fID = open(string(folder,"/result.dat"),"w")
    println(fID,string("T: ",mean(T)," K"))
    println(fID,string("ρ_v: ",ρ_v," g/ml"))
    println(fID,string("ρ_l: ",ρ_l," g/ml"))
    println(fID,string("p_v: ",p_v," MPa"))
    println(fID,string("p_l: ",p_l," MPa"))
    println(fID,string("D: ",D," Å"))
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

# Function to read profile data
function read_profile1D(filename,data,ts_add)
    fID = open(filename,"r")

    # Read first three lines
    line1 = readline(fID)
    line2 = readline(fID)
    line3 = readline(fID)
    if (line2 != "# Timestep Number-of-chunks Total-count") ||
       (line3 != "# Chunk Coord1 Ncount density/number density/mass temp v_p_xx v_p_yy v_p_zz v_p_xy v_p_xz v_p_yz")
        error("Dataset not fitting to 'load_profile1D' function!")
    end

    # Read data body
    txt = readlines(fID)

    # Determine number of chunks
    pos = vcat(findall(isequal(' '),txt[1]),length(txt[1])+1)
    no_chunks = parse(Int64,txt[1][pos[1]+1:pos[2]-1])
    lines_per_ts = (no_chunks+1)
    n_steps = Int64(length(txt)/lines_per_ts)

    if typeof(data) != profile_data
        # Initialization of profile_data strucutre
        init = Array{Float64,2}(undef,0,no_chunks)
        data = profile_data(Float64[],init,init,init,init,init,init,init,init,init,init,init,init)
    end

    for i = 1:n_steps
        # Get lines of step
        istart = lines_per_ts*(i-1) + 1
        iend = istart + no_chunks

        # Get information from first line
        line_start = txt[istart]
        pos = vcat(findall(isequal(' '),line_start),length(line_start)+1)
        timestep = parse(Int64,line_start[1:pos[1]-1]) + ts_add
        total_count = parse(Float64,line_start[pos[2]+1:end])
        append!(data.timestep,timestep)

        # cols = length(properties)
        body = Array(transpose(parse.(Float64,hcat(split.(txt[istart+1:iend])...))))
        data.id_chunk = vcat(data.id_chunk,body[:,1]')
        data.x = vcat(data.x,body[:,2]')
        data.Ncount = vcat(data.Ncount,body[:,3]')
        data.ρn = vcat(data.ρn,body[:,4]')
        data.ρm = vcat(data.ρm,body[:,5]')
        data.T = vcat(data.T,body[:,6]')
        data.pxx = vcat(data.pxx,body[:,7]').*0.1
        data.pyy = vcat(data.pyy,body[:,8]').*0.1
        data.pzz = vcat(data.pzz,body[:,9]').*0.1
        data.pxy = vcat(data.pxy,body[:,10]').*0.1
        data.pxz = vcat(data.pxz,body[:,11]').*0.1
        data.pyz = vcat(data.pyz,body[:,12]').*0.1
    end

    return data
end

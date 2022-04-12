## LoadData.jl
# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - LoadData
# Function to load data from simulation folder to variables
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

# Functions to load different files
# Loading Info File
function load_info(folder; is_nemd="no")
    # Get filename
    list = readdir(folder)
    filename = list[occursin.("info.",list)][end]
    file = string(folder,"/",filename)

    if is_nemd == "no"
        # Read file info.dat
        if isfile(file)
            fID = open(file,"r");   lines = readlines(fID);     close(fID)

            # Extract information
            pos1 = findfirst(": ",lines[1])
            moltype = lines[1][pos1[end]+1:end]
            pos2 = findfirst(": ",lines[2])
            dt = parse(Float64,lines[2][pos2[end]+1:end])
            pos3 = findfirst(": ",lines[3])
            natoms = parse(Int64,lines[3][pos3[end]+1:end])
            pos4 = findfirst(": ",lines[4])
            molmass = parse(Float64,lines[4][pos4[end]+1:end])

            return moltype, dt, natoms, molmass

        else error("File \"",file,"\" is empty") end
    elseif is_nemd == "shear"
    elseif is_nemd == "heat"
        # Read file info.dat
        if isfile(file)
            fID = open(file,"r");   lines = readlines(fID);     close(fID)

            # Extract information
            pos5 = findfirst(": ",lines[5])
            Lx = parse(Float64,lines[5][pos5[end]+1:end])
            pos6 = findfirst(": ",lines[6])
            Ly = parse(Float64,lines[6][pos6[end]+1:end])
            pos7 = findfirst(": ",lines[7])
            Lz = parse(Float64,lines[7][pos7[end]+1:end])
            pos8 = findfirst(": ",lines[8])
            nbins = parse(Float64,lines[8][pos8[end]+1:end])

            return Lx, Ly, Lz

        else error("File \"",file,"\" is empty") end
    end
end

# Loading Thermo File
function load_thermo(info::info_struct; is_nemd::String="no")
    thermodat = thermo_dat(Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[])
    list = sort(readdir(info.folder))
    files = string.(info.folder,"/",["LJTrans_out.res"])
    # if length(files) > 9
    #     files = files[sortperm(parse.(Int64,getindex.(split.(files,"."),2)))]
    # end
    stepADD=0
    timeADD=0
    for file in files
        if is_nemd == "no"
            step, time, T, p, ρ, Etot, Ekin, Epot = load_thermo_file(file,info,is_nemd)
        elseif is_nemd == "shear"
            step, time, T, p, ρ, Etot, Ekin, Epot, pyz = load_thermo_file(file,info,is_nemd)
        elseif is_nemd == "heat"
            step, time, T, p, ρ, Etot, Ekin, Epot, Qhot, Qcold = load_thermo_file(file,info,is_nemd)
        end
        thermodat.step = vcat(thermodat.step,step.+stepADD)
        thermodat.t = vcat(thermodat.t,time.+timeADD)
        thermodat.T = vcat(thermodat.T,T)
        thermodat.p = vcat(thermodat.p,p)
        thermodat.ρ = vcat(thermodat.ρ,ρ)
        thermodat.Etot = vcat(thermodat.Etot,Etot)
        thermodat.Ekin = vcat(thermodat.Ekin,Ekin)
        thermodat.Epot = vcat(thermodat.Epot,Epot)
        if is_nemd == "shear"
            thermodat.pyz = vcat(thermodat.pyz,pyz)
        elseif is_nemd == "heat"
            thermodat.Qhot = vcat(thermodat.Qhot,Qhot)
            thermodat.Qcold = vcat(thermodat.Qcold,Qcold)
        end
        stepADD = thermodat.step[end]
        timeADD = thermodat.t[end]
    end

    if isempty(thermodat.step) error("No thermo data loaded!") end
    return thermodat
end

function load_thermo_file(file::String, info::info_struct, is_nemd::String)
    fID = open(file,"r"); readline(fID); readline(fID); readline(fID); line4 = readline(fID); close(fID)

    if is_nemd == "no"
        if line4 != "#step	t		U_pot	U_pot_avg		p	p_avg		beta_trans	beta_rot		c_v		N"
            error("Format of File \"LJTrans_out.res\" not right")
        end
        # Read data
        dat = readdlm(file, skipstart=4)
        step = Int64.(dat[1:end-1,1]);  time = Float64.(dat[1:end-1,2])
        T    = Float64[];               p    = Float64.(dat[1:end-1,5])
        ρ    = Float64[];               Etot = Float64[]
        Ekin = Float64[];               Epot = Float64.(dat[1:end-1,3])
        return step, time, T, p, ρ, Etot, Ekin, Epot

    elseif is_nemd == "shear"
        if line2 != "# TimeStep v_T v_p v_rho v_Etot v_Ekin v_Epot v_pyz v_eta"
            error("Format of File \"thermo.dat\" not right")
        end
        # Read data
        dat = readdlm(file, skipstart=2)
        step = dat[:,1];    time = step*info.dt
        T    = dat[:,2];    p    = dat[:,3]
        ρ    = dat[:,4];    Etot = dat[:,5]
        Ekin = dat[:,6];    Epot = dat[:,7]
        pyz  = dat[:,8]
        return step, time, T, p, ρ, Etot, Ekin, Epot, pyz

    elseif is_nemd == "heat"
        if line2 != "# TimeStep v_T v_p v_rho v_Etot v_Ekin v_Epot f_hot f_cold"
            error("Format of File \"thermo.dat\" not right")
        end
        # Read data
        dat = readdlm(file, skipstart=2)
        step = dat[:,1];    time = step*info.dt
        T    = dat[:,2];    p    = dat[:,3]
        ρ    = dat[:,4];    Etot = dat[:,5]
        Ekin = dat[:,6];    Epot = dat[:,7]
        Qhot = dat[:,8];    Qcold= dat[:,9]
        return step, time, T, p, ρ, Etot, Ekin, Epot, Qhot, Qcold
    end
end

# Loading Pressure File
function load_pressure(info, nEvery)
    pdat = pressure_dat(Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[])
    list = sort(readdir(info.folder))
    files = string.(info.folder,"/",["PressureTensor.Vipr"])
    # if length(files) > 9
    #     files = files[sortperm(parse.(Int64,getindex.(split.(files,"."),2)))]
    # end
    let stepADD=0,timeADD=0, skip=0
        for file in files
            step, time, V, T, p, pxx, pyy, pzz, pxy, pxz, pyz = load_pressure_file(file, info, skip, nEvery)
            pdat.step = vcat(pdat.step,step.+stepADD)
            pdat.t = vcat(pdat.t,time.+timeADD)
            pdat.V = vcat(pdat.V,V)
            pdat.T = vcat(pdat.T,T)
            pdat.p = vcat(pdat.p,p)
            pdat.pxx = vcat(pdat.pxx,pxx)
            pdat.pyy = vcat(pdat.pyy,pyy)
            pdat.pzz = vcat(pdat.pzz,pzz)
            pdat.pxy = vcat(pdat.pxy,pxy)
            pdat.pxz = vcat(pdat.pxz,pxz)
            pdat.pyz = vcat(pdat.pyz,pyz)
            stepADD = pdat.step[end]
            timeADD = pdat.t[end]
            if (skip == 0) skip = 1 end
        end
    end
    if isempty(pdat.step)
        pdat = []
        @warn("No pressure data loaded!")
    end
    return pdat
end

function load_pressure_file(file, info, skip, nEvery)
    lines = readlines(file)
    if lines[12] != "# old prefix string (get time step here)	globalT	p_xx^[1]	p_yy^[1]	p_zz^[1]	p_xx^[2]	p_yy^[2]	p_zz^[2]	p_xy^[2]	p_xz^[2]	p_yz^[2]"
        error("Format of File \"PressureTensor.Vipr\" not right")
    end
    Lxyz = parse.(Float64,split(lines[5]))
    V = Lxyz[1]*Lxyz[2]*Lxyz[3]
    # Read data
    dat = readdlm(file, skipstart=(13+skip))
    what = skip*nEvery+1-skip:nEvery:size(dat,1)
    if mod(size(dat,1)-1+skip,nEvery) != 0
        error("n_every must be a divisor of $(size(dat,1)-1+skip)!")
    end
    step = str2timestep.(dat[what,1]);      step = step .- step[1]
    time = step*info.dt
    V    = V.*ones(length(step));    T    = dat[what,2];    p    = Float64[]
    pxx  = dat[what,6];    pyy  = dat[what,7];    pzz  = dat[what,8]
    pxy  = dat[what,9];    pxz  = dat[what,10];    pyz  = dat[what,11]
    return step, time, V, T, p, pxx, pyy, pzz, pxy, pxz, pyz
end

function str2timestep(in)
    out = parse(Int64,split(in,"_")[2])
end

# Loading Dump File
function load_dump(info)
    # Initialization of posdat
    posdat = []

    # Get list of all atoms_position files
    list = sort(readdir(info.folder))
    files = string.(info.folder,"/",list[startswith.(list,"LJTrans_xyz_-")])
    ts_file = get_timestep_from_file.(files).*5000
    sort_id = sortperm(ts_file)
    files = files[sort_id]
    ts_file = ts_file[sort_id]

    # Read dump file timestep by timestep
    let stepADD=0, timeADD=0, skip=0, i = 0
        for file in files
            i += 1
            fID = open(file,"r");   readline(fID);  line2 = readline(fID); close(fID)
            if line2 == "comment line"
                dat  = readdlm(file, skipstart=2)
                if !(ts_file[i]==0 && skip==1)
                    posdat = vcat(posdat,dump_dat(ts_file[i], ts_file[i]*info.dt, [0.0 1.0;0.0 1.0;0.0 1.0], dat[:,1], ones(Int64,info.natoms), dat[:,1], ones(Int64,info.natoms), ones(info.natoms), dat[:,2], dat[:,3], dat[:,4]))
                end
            else
                error("Wrong file format for ls1!")
            end
            skip = 1
            stepADD = posdat[end].step
            timeADD = posdat[end].t
        end
    end

    # Sort data such that ID of atoms are sorted
    for i = 1:length(posdat)
        p = posdat[i]
        i_sort = sortperm(p.id)
        for f in fieldnames(typeof(p))
            if !(f in [:step, :t, :bounds, :moltype])
                setfield!(p,f,getfield(p,f)[i_sort])
            end
        end
        posdat[i] = p
    end

    return posdat
end

function get_timestep_from_file(filename)
    ts = parse(Int64,split(filename,"xyz")[2][3:end-1])
end

# Loading Heat Flux File
function load_heatflux(info, nEvery)
    jdat = heat_dat(Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[])
    list = sort(readdir(info.folder))
    files = string.(info.folder,"/",list[occursin.(string("heat_flux.",info.ensemble,"."),list)])
    if length(files) > 9
        files = files[sortperm(parse.(Int64,getindex.(split.(files,"."),2)))]
    end
    let stepADD=0,timeADD=0, skip=0
        for file in files
            step, time, V, T, jx, jy, jz = load_heatflux_file(file, info, skip, nEvery)
            jdat.step = vcat(jdat.step,step.+stepADD)
            jdat.t = vcat(jdat.t,time.+timeADD)
            jdat.V = vcat(jdat.V,V)
            jdat.T = vcat(jdat.T,T)
            jdat.jx = vcat(jdat.jx,jx)
            jdat.jy = vcat(jdat.jy,jy)
            jdat.jz = vcat(jdat.jz,jz)
            stepADD = jdat.step[end]
            timeADD = jdat.t[end]
            if (skip == 0) skip = 1 end
        end
    end
    if isempty(jdat.step)
        jdat = []
        @warn("No heatflux data loaded!")
    end
    return jdat
end
function load_heatflux_file(file, info, skip, nEvery)
    fID = open(file,"r"); readline(fID); line2 = readline(fID); close(fID)
    if line2 != "# TimeStep v_V v_T c_j[1] c_j[2] c_j[3]"
        error("Format of File \"heat_flux.dat\" not right")
    end
    # Read data
    dat = readdlm(file, skipstart=(2+skip))
    what = skip*nEvery+1-skip:nEvery:size(dat,1)
    if mod(size(dat,1)-1+skip,nEvery) != 0
        error("n_every must be a divisor of $(size(dat,1)-1+skip)!")
    end
    step = dat[what,1];    time = step*info.dt
    V    = dat[what,2];    T    = dat[what,3]
    jx   = dat[what,4];    jy   = dat[what,5];    jz = dat[what,6]
    return step, time, V, T, jx, jy, jz
end

# Loading thermo file in vle simulations
function load_thermo_vle(folder::String, info::info_struct)
    thermodat = thermo_vle_dat(Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[])
    list = sort(readdir(folder))
    files = string.(folder,"/",list[occursin.(string("thermo.vle.2phase."),list)])
    if length(files) > 9
        files = files[sortperm(parse.(Int64,getindex.(split.(files,"."),2)))]
    end
    stepADD=0
    timeADD=0
    for file in files
        step, time, T, px, py, pz, ρ, Etot, Ekin, Epot = load_thermo_vle_file(file,info)
        thermodat.step = vcat(thermodat.step,step.+stepADD)
        thermodat.t = vcat(thermodat.t,time.+timeADD)
        thermodat.T = vcat(thermodat.T,T)
        thermodat.px = vcat(thermodat.px,px)
        thermodat.py = vcat(thermodat.py,py)
        thermodat.pz = vcat(thermodat.pz,pz)
        thermodat.ρ = vcat(thermodat.ρ,ρ)
        thermodat.Etot = vcat(thermodat.Etot,Etot)
        thermodat.Ekin = vcat(thermodat.Ekin,Ekin)
        thermodat.Epot = vcat(thermodat.Epot,Epot)
        stepADD = thermodat.step[end]
        timeADD = thermodat.t[end]
    end

    return thermodat
    if isempty(thermodat.step) error("No thermo data loaded!") end
    return thermodat
end
function load_thermo_vle_file(file::String, info::info_struct)
    fID = open(file,"r"); readline(fID); line2 = readline(fID); close(fID)
    if line2 != "# TimeStep v_T v_px v_py v_pz v_rho v_Etot v_Ekin v_Epot"
        error("Format of File \"thermo.vle.2phase.dat\" not right")
    end
    # Read data
    dat = readdlm(file, skipstart=2)
    step = dat[:,1];    time = step*info.dt
    T    = dat[:,2];    px    = dat[:,3]
    py    = dat[:,4];   pz    = dat[:,5]
    ρ    = dat[:,6];    Etot = dat[:,7]
    Ekin = dat[:,8];    Epot = dat[:,9]
    return step, time, T, px, py, pz, ρ, Etot, Ekin, Epot
end

# Function to read profile data
function read_profile1D(filename,data,ts_add)
    fID = open(filename,"r")

    # Read first three lines
    line1 = readline(fID)
    line2 = readline(fID)
    line3 = readline(fID)
    type = ""
    if line2 == "# Timestep Number-of-chunks Total-count"
        if line3 == "# Chunk Coord1 Ncount density/number density/mass temp v_p_xx v_p_yy v_p_zz v_p_xy v_p_xz v_p_yz"
            type = "vle"
        elseif line3 == "# Chunk Coord1 Ncount vy"
            type = "shear"
        elseif line3 == "# Chunk Coord1 Ncount temp"
            type = "heat"
        end
    elseif line2 == "# TimeStep Number-of-rows"
        if startswith(line3,"# Row c_myRDF[1] c_myRDF[2] c_myRDF[3]")
            N_rdf = Int64((sum(startswith.(split(line3),"c_myRDF")) - 1)/2)
            type = "rdf"
        end
    end
    if isempty(type) error("File format of \"$filename\" is unknown to 'load_profile1D' function!") end

    # Read data body
    txt = readlines(fID)

    # Determine number of chunks
    pos = vcat(findall(isequal(' '),txt[1]),length(txt[1])+1)
    no_chunks = parse(Int64,txt[1][pos[1]+1:pos[2]-1])
    lines_per_ts = (no_chunks+1)
    n_steps = Int64(length(txt)/lines_per_ts)

    if (typeof(data) != profile_data_vle) || (typeof(data) != profile_data_shear) || (typeof(data) != profile_data_heat) || (typeof(data) != profile_data_rdf)
        # Initialization of profile_data strucutre
        init = Array{Float64,2}(undef,0,no_chunks)
        if type == "vle"
            data = profile_data_vle(Float64[],init,init,init,init,init,init,init,init,init,init,init,init)
        elseif type == "shear"
            data = profile_data_shear(Float64[],init,init,init,init)
        elseif type == "heat"
            data = profile_data_heat(Float64[],init,init,init,init)
        elseif type == "rdf"
            data = profile_data_rdf(Float64[],init,init,repeat([init],N_rdf),repeat([init],N_rdf))
        end
    end

    for i = 1:n_steps
        # Get lines of step
        istart = lines_per_ts*(i-1) + 1
        iend = istart + no_chunks

        # Get information from first line
        line_start = txt[istart]
        timestep = parse(Int64,split(line_start)[1]) + ts_add
        append!(data.timestep,timestep)

        # cols = length(properties)
        body = Array(transpose(parse.(Float64,hcat(split.(txt[istart+1:iend])...))))
        data.id_chunk = vcat(data.id_chunk,body[:,1]')
        if type == "vle"
            data.x = vcat(data.x,body[:,2]')
            data.Ncount = vcat(data.Ncount,body[:,3]')
            data.ρn = vcat(data.ρn,body[:,4]')
            data.ρm = vcat(data.ρm,body[:,5]')
            data.T = vcat(data.T,body[:,6]')
            data.pxx = vcat(data.pxx,body[:,7]')
            data.pyy = vcat(data.pyy,body[:,8]')
            data.pzz = vcat(data.pzz,body[:,9]')
            data.pxy = vcat(data.pxy,body[:,10]')
            data.pxz = vcat(data.pxz,body[:,11]')
            data.pyz = vcat(data.pyz,body[:,12]')
        elseif type == "shear"
            data.x = vcat(data.x,body[:,2]')
            data.Ncount = vcat(data.Ncount,body[:,3]')
            data.vy = vcat(data.vy,body[:,4]')
        elseif type == "heat"
            data.x = vcat(data.x,body[:,2]')
            data.Ncount = vcat(data.Ncount,body[:,3]')
            data.T = vcat(data.T,body[:,4]')
        elseif type == "rdf"
            data.r = vcat(data.r,body[:,2]')
            for i in 1:N_rdf
                data.g_r[i] = vcat(data.g_r[i],body[:,3+(i-1)*2]')
                data.N_coord[i] = vcat(data.N_coord[i],body[:,4+(i-1)*2]')
            end
        end
    end

    return data, no_chunks, n_steps
end

# Function to read molcount/atoms per molecule
function read_inputfile(filename)
    fID = open(filename,"r")

    # Read first three lines
    #line1 = readline(fID)
    #line2 = readline(fID)
    #line3 = readline(fID)
    #type = ""
    #if line2 == "variable        N           equal"
    txts=readlines(fID)
    tst=length(txts)
    close(fID)
    fID = open(filename,"r")
    local N_mol, n
    for i in 1:tst
        txt=readline(fID)
        string_0 = replace(txt," " => "" )
        if startswith(string_0, "variableNequal")
            n=findfirst("variableNequal", string_0)
            N_mol=string_0[15:end]
            m=findfirst("#", N_mol)
            N_mol=N_mol[1:m[1]-1]
            N_mol=parse(Int,N_mol)
        end
    end
    return N_mol
end

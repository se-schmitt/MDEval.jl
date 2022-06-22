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
    files = string.(info.folder,"/",list[occursin.(string("thermo.",info.ensemble,"."),list)])
    if length(files) > 9
        files = files[sortperm(parse.(Int64,getindex.(split.(files,"."),2)))]
    end
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
    if is_nemd == "gro"
    else
    fID = open(file,"r"); readline(fID); line2 = readline(fID); close(fID)
    end
    if is_nemd == "no"
        if line2 != "# TimeStep v_T v_p v_rho v_Etot v_Ekin v_Epot"
            error("Format of File \"thermo.dat\" not right")
        end
        # Read data
        dat = readdlm(file, skipstart=2)
        step = dat[:,1];    time = step*info.dt
        T    = dat[:,2];    p    = dat[:,3]
        ρ    = dat[:,4];    Etot = dat[:,5]
        Ekin = dat[:,6];    Epot = dat[:,7]
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
    pdat = pressure_gro(Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[])
    list = sort(readdir(info.folder))
    files = string.(info.folder,"/",list[occursin.(string("pressure.",info.ensemble,"."),list)])
    if length(files) > 9
        files = files[sortperm(parse.(Int64,getindex.(split.(files,"."),2)))]
    end
    file = string.(info.folder,"/","density.xvg")
    #file = "C:/Users/kn9-f/md-evaluation/T0.5rho0.9s0.1/density.xvg"
    #file = "C:/Users/kn9-f/md-evaluation/T1.56rho0.546s0.1/density.xvg"
    #file = "C:/Users/kn9-f/md-evaluation/T1.56rho0.1075s0.1/density.xvg"
    let stepADD=0,timeADD=0, skip=0
            time, Epot, Ekin, Etot, T, p, pxx, pxy, pxz, pyy, pyz, pzz = load_pressure_file(file, info, skip, nEvery)
            pdat.t = vcat(pdat.t,time.+timeADD)
            pdat.Epot = vcat(pdat.Epot,Epot)
            pdat.Ekin = vcat(pdat.Ekin,Ekin)
            pdat.Etot = vcat(pdat.Etot,Etot)
            pdat.T = vcat(pdat.T,T)
            pdat.p = vcat(pdat.p,p)
            pdat.pxx = vcat(pdat.pxx,pxx)
            pdat.pxy = vcat(pdat.pxy,pxy)
            pdat.pxz = vcat(pdat.pxz,pxz)
            pdat.pyy = vcat(pdat.pyy,pyy)
            pdat.pyz = vcat(pdat.pyz,pyz)
            pdat.pzz = vcat(pdat.pzz,pzz)
            #stepADD = pdat.step[end]
            #timeADD = pdat.t[end]
            if (skip == 0) skip = 1 end
    end
    if isempty(pdat.t)
        pdat = []
        @warn("No pressure data loaded!")
    end
    return pdat
end

function load_pressure_file(file, info, skip, nEvery)
    fID = open(file,"r"); line=readline(fID); skip_gro=1; end_gro=1;
    while end_gro == 1
        if startswith(line, "#")
            line = readline(fID);
            skip_gro = skip_gro + 1;
        elseif startswith(line, "@")
            line = readline(fID);
            skip_gro = skip_gro + 1;
        else
            end_gro = 0;
            skip_gro = skip_gro - 1;
        end
    end
    close(fID)
    # Read data
    dat = readdlm(file, skipstart=(skip_gro+skip))
    what = skip*nEvery+1-skip:nEvery:size(dat,1)
    if mod(size(dat,1)-1+skip,nEvery) != 0
        error("n_every must be a divisor of $(size(dat,1)-1+skip)!")
    end
    fac_gro_p = 1 ./ 16.60544783
    #fac_gro_p = 1
    time = dat[what,1];
    Epot = dat[what,2];    Ekin = dat[what,3];    Etot = dat[what,4]
    T    = dat[what,6] ./120.27236;
    p    = dat[what,7] .* fac_gro_p;    pxx  = dat[what,8] .* fac_gro_p;    pxy  = dat[what,9] .* fac_gro_p
    pxz  = dat[what,10] .* fac_gro_p;   pyy  = dat[what,11] .* fac_gro_p;   pyz  = dat[what,12] .* fac_gro_p
    pzz  = dat[what,13] .* fac_gro_p
    return time, Epot, Ekin, Etot, T, p, pxx, pxy, pxz, pyy, pyz, pzz
end

# Loading Dump File
function load_dump(info)
    # Initialization of posdat
    posdat = []

    # Get list of all atoms_position files
    list = sort(readdir(info.folder))
    files = string.(info.folder,"/",list[occursin.(string("atoms_position.",info.ensemble,"."),list)])
    if length(files) > 9
        files = files[sortperm(parse.(Int64,getindex.(split.(files,"."),2)))]
    end

    # Read dump file timestep by timestep
    let stepADD=0, timeADD=0, skip=0
        for file in files
            fID = open(file,"r")
            while !eof(fID)
                # Reading timestep and calculating time
                if (readline(fID) == "ITEM: TIMESTEP")
                    step = parse(Int64,readline(fID))
                    t = step*info.dt
                else error("Dump file format wrong") end
                # Reading number of atoms
                if (readline(fID) == "ITEM: NUMBER OF ATOMS")
                    natoms = parse(Int16,readline(fID))
                else error("Dump file format wrong") end
                # Reading boundary coords
                if (readline(fID) == "ITEM: BOX BOUNDS pp pp pp")
                    bounds = zeros(3,2)
                    bounds[1,:] = parse.(Float64,split(readline(fID)))
                    bounds[2,:] = parse.(Float64,split(readline(fID)))
                    bounds[3,:] = parse.(Float64,split(readline(fID)))
                else error("Dump file format wrong") end
                # Reading positions of atoms
                lineITEM = readline(fID)
                id = Int64.(zeros(natoms))
                molid = Int64.(zeros(natoms))
                type = Int64.(zeros(natoms))
                x = zeros(natoms)
                y = zeros(natoms)
                z = zeros(natoms)
                if (lineITEM == "ITEM: ATOMS id mol xu yu zu ")
                    mass = []
                    for i = 1:natoms
                        line_float = parse.(Float64,split(readline(fID)))
                        id[i] = Int64(line_float[1])
                        molid[i] = Int64(line_float[2])
                        x[i] = line_float[3]
                        y[i] = line_float[4]
                        z[i] = line_float[5]
                    end
                elseif (lineITEM == "ITEM: ATOMS id mol mass xu yu zu ")
                    mass = zeros(natoms)
                    for i = 1:natoms
                        line_float = parse.(Float64,split(readline(fID)))
                        id[i] = Int64(line_float[1])
                        molid[i] = Int64(line_float[2])
                        mass[i] = line_float[3]
                        x[i] = line_float[4]
                        y[i] = line_float[5]
                        z[i] = line_float[6]
                    end
                elseif (lineITEM == "ITEM: ATOMS id mass xu yu zu ")
                    mass = zeros(natoms)
                    for i = 1:natoms
                        line_float = parse.(Float64,split(readline(fID)))
                        id[i] = Int64(line_float[1])
                        molid[i] = id[i]
                        mass[i] = line_float[2]
                        x[i] = line_float[3]
                        y[i] = line_float[4]
                        z[i] = line_float[5]
                    end
                elseif (lineITEM == "ITEM: ATOMS id mol type mass xu yu zu ")
                    mass = zeros(natoms)
                    for i = 1:natoms
                        line_float = parse.(Float64,split(readline(fID)))
                        id[i] = Int64(line_float[1])
                        molid[i] = Int64(line_float[2])
                        type[i] = Int64(line_float[3])
                        mass[i] = line_float[4]
                        x[i] = line_float[5]
                        y[i] = line_float[6]
                        z[i] = line_float[7]
                    end
                else error("Dump file format wrong") end

                if !(step==0 && skip==1)
                    posdat = vcat(posdat,dump_dat(step+stepADD, t+timeADD, bounds, id, type, molid, Int64[], mass, x, y, z))
                end
            end
            skip = 1
            stepADD = posdat[end].step
            timeADD = posdat[end].t
        end
    end

    if !(isempty(posdat))
        # Caluclation of moltypes (just for first timestep)
        Nmol = length(unique(posdat[1].molid))
        moltype = zeros(Int64,Nmol)
        types_mol = []

        # Loop over all molecules
        i = 0
        for molid in sort(unique(posdat[1].molid))
            i += 1
            what = posdat[1].molid .== molid
            types = sort(posdat[1].type[what])

            if types in types_mol
                moltype[i] = findfirst(types_mol .== [types])
            else
                moltype[i] = maximum(moltype)+1
                append!(types_mol,[types])
            end
        end
        posdat[1].moltype = moltype
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

function load_gro_file(file::String, info::info_struct)
    if file == "MSD.xvg" # !! "MSD.xvg"
    file = string(info.folder,"/",file)
    msd = String[]
    fID = open(file,"r"); readline(fID); readline(fID); line3 = readline(fID); close(fID)
    if line3 != "#                        :-) GROMACS - gmx msd, 2022 (-:"
        error("Format of File \"MSD.xvg\" not right")
    end
    for line in eachline(file)
        if startswith(line,"@ s0 legend")
            msd = line
        end
    end
    msd1 = parse(Float64,SubString(msd,(findfirst("=",msd))[1]+2,(findfirst("(",msd))[1]-2))
    msd2 = parse(Float64,SubString(msd,(findfirst("-",msd))[1]+2,(findfirst(")",msd))[1]-1))
    return msd1, msd2
    end
end

function load_nemd_file(info::info_struct)
    list = sort(readdir(info.folder))
    nemd_string = "density_nemd_"
    files = string.(info.folder,"/",list[occursin.(string("density_nemd_"),list)])
    s_rates = zeros(length(files))
    for i = 1:length(files)
        substring = files[i]
        a=findlast(nemd_string,substring)
        b=findlast(".xvg",substring)
        a=a[end]
        b=b[1]
        s_rates[i]=parse(Float64, substring[a+1:b-1])
    end
    s_rate = maximum(s_rates)
    filename = string(nemd_string,s_rate,".xvg")

    file = string(info.folder,"/",filename)
    fID = open(file,"r"); line=readline(fID); skip_gro=1; end_gro=1;
    while end_gro == 1
        if startswith(line, "#")
            line = readline(fID);
            skip_gro = skip_gro + 1;
        elseif startswith(line, "@")
            line = readline(fID);
            skip_gro = skip_gro + 1;
        else
            end_gro = 0;
            skip_gro = skip_gro - 1;
        end
    end
    close(fID)

    # Read data
    dat = readdlm(file, skipstart=skip_gro)
    what = 1:1:size(dat,1)

    time = dat[what,1];
    vis0 = dat[what,2];
    return time, vis0, s_rate
end

function load_etta_files(info::info_struct)
    #load gromacs files for Einstein viscosity calc.
    list = sort(readdir(info.folder))
    etta_string = "evisco"
    files = string.(info.folder,"/",list[occursin.(string("evisco"),list)])
    counter = 1
    let η_t_all=0.0, t=0.0, counter = 1
    for i=1:length(files)
        file=files[i]
        #file = string(info.folder,"/",filename)
        fID = open(file,"r"); line=readline(fID); skip_gro=1; end_gro=1;
        while end_gro == 1
            if startswith(line, "#")
                line = readline(fID);
                skip_gro = skip_gro + 1;
            elseif startswith(line, "@")
                line = readline(fID);
                skip_gro = skip_gro + 1;
            else
                end_gro = 0;
                skip_gro = skip_gro - 1;
            end
        end
        close(fID)

    # Read data
    dat = readdlm(file, skipstart=skip_gro)
    η_t_app = dat[1:end,2:end]
    if counter == 1
        t=dat[1:end,1]
        η_t_all = η_t_app
        counter = 0
    else
        η_t_all = hcat(η_t_all,η_t_app)
    end
    end
    return t, η_t_all
    end

    @infiltrate
    return t, η_t_all
end

# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - LoadData
# Function to load data from simulation folder to variables
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

# Functions to load different files
# Loading Info File
function load_info(folder)
    # Get filename
    list = readdir(folder)
    filename = list[occursin.("info.",list)][end]
    file = string(folder,"/",filename)

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
    else error("File \"",file,"\" is empty") end

    return moltype, dt, natoms, molmass
end

# Loading Thermo File
function load_thermo(info::info_struct)
    thermodat = thermo_dat(Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[])
    list = sort(readdir(info.folder))
    files = string.(info.folder,"/",list[occursin.(string("thermo.",info.ensemble,"."),list)])
    if length(files) > 9
        files = files[sortperm(parse.(Int64,getindex.(split.(files,"."),2)))]
    end
    stepADD=0
    timeADD=0
    for file in files
        step, time, T, p, ρ, Etot, Ekin, Epot = load_thermo_file(file,info)
        thermodat.step = vcat(thermodat.step,step.+stepADD)
        thermodat.t = vcat(thermodat.t,time.+timeADD)
        thermodat.T = vcat(thermodat.T,T)
        thermodat.p = vcat(thermodat.p,p)
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
function load_thermo_file(file::String, info::info_struct)
    fID = open(file,"r"); readline(fID); line2 = readline(fID); close(fID)
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
end

# Loading Pressure File
function load_pressure(info)
    pdat = pressure_dat(Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[])
    list = sort(readdir(info.folder))
    files = string.(info.folder,"/",list[occursin.(string("pressure.",info.ensemble,"."),list)])
    if length(files) > 9
        files = files[sortperm(parse.(Int64,getindex.(split.(files,"."),2)))]
    end
    let stepADD=0,timeADD=0, skip1=0
        for file in files
            step, time, V, T, p, pxx, pyy, pzz, pxy, pxz, pyz = load_pressure_file(file,info,skip1)
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
            if (skip1 == 0) skip1 = 1 end
        end
    end
    if isempty(pdat.step)
        pdat = []
        @warn("No pressure data loaded!")
    end
    return pdat
end

function load_pressure_file(file,info,skip)
    fID = open(file,"r"); readline(fID); line2 = readline(fID); close(fID)
    if line2 != "# TimeStep v_V v_T c_thermo_press c_thermo_press[1] c_thermo_press[2] c_thermo_press[3] c_thermo_press[4] c_thermo_press[5] c_thermo_press[6]"
        error("Format of File \"pressure.dat\" not right")
    end
    # Read data
    dat = readdlm(file, skipstart=(2+skip))
    step = dat[:,1];    time = step*info.dt
    V    = dat[:,2];    T    = dat[:,3];    p    = dat[:,4]
    pxx  = dat[:,5];    pyy  = dat[:,6];    pzz  = dat[:,7]
    pxy  = dat[:,8];    pxz  = dat[:,9];    pyz  = dat[:,10]
    return step, time, V, T, p, pxx, pyy, pzz, pxy, pxz, pyz
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
    let stepADD=0, timeADD=0, skip1=0
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

                if !(step==0 && skip1==1)
                    posdat = vcat(posdat,dump_dat(step+stepADD, t+timeADD, bounds, id, type, molid, Int64[], mass, x, y, z))
                end
            end
            skip1 = 1
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

    return posdat
end

# Loading Heat Flux File
function load_heatflux(info)
    jdat = heat_dat(Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[])
    list = sort(readdir(info.folder))
    files = string.(info.folder,"/",list[occursin.(string("heat_flux.",info.ensemble,"."),list)])
    if length(files) > 9
        files = files[sortperm(parse.(Int64,getindex.(split.(files,"."),2)))]
    end
    let stepADD=0,timeADD=0, skip1=0
        for file in files
            step, time, V, T, jx, jy, jz = load_heatflux_file(file,info,skip1)
            jdat.step = vcat(jdat.step,step.+stepADD)
            jdat.t = vcat(jdat.t,time.+timeADD)
            jdat.V = vcat(jdat.V,V)
            jdat.T = vcat(jdat.T,T)
            jdat.jx = vcat(jdat.jx,jx)
            jdat.jy = vcat(jdat.jy,jy)
            jdat.jz = vcat(jdat.jz,jz)
            stepADD = jdat.step[end]
            timeADD = jdat.t[end]
            if (skip1 == 0) skip1 = 1 end
        end
    end
    if isempty(jdat.step)
        jdat = []
        @warn("No heatflux data loaded!")
    end
    return jdat
end
function load_heatflux_file(file,info,skip)
    fID = open(file,"r"); readline(fID); line2 = readline(fID); close(fID)
    if line2 != "# TimeStep v_V v_T c_j[1] c_j[2] c_j[3]"
        error("Format of File \"heat_flux.dat\" not right")
    end
    # Read data
    dat = readdlm(file, skipstart=(2+skip))
    step = dat[:,1];    time = step*info.dt
    V    = dat[:,2];    T    = dat[:,3]
    jx   = dat[:,4];    jy   = dat[:,5];    jz = dat[:,6]
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
        data.pxx = vcat(data.pxx,body[:,7]')
        data.pyy = vcat(data.pyy,body[:,8]')
        data.pzz = vcat(data.pzz,body[:,9]')
        data.pxy = vcat(data.pxy,body[:,10]')
        data.pxz = vcat(data.pxz,body[:,11]')
        data.pyz = vcat(data.pyz,body[:,12]')
    end

    return data
end

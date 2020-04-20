# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - LoadData
# Function to load data from simulation folder to variables
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

# Main Function
function LoadData(info)

    # Loading Info File
    info = load_info(info)

    # Loading Thermo File
    thermodat = load_thermo(info)

    # Loading Pressure File
    pdat = load_pressure(info)

    # Loading Dump File
    # posdat = load_dump(info)
    posdat = []

    # Loading Heat Flux File
    # jdat = load_heatflux(info)
    jdat = []

    DATA = dats2DATA(thermodat, pdat, posdat, jdat)
    return DATA, info
end

# Subfunctions
# Loading Info File
function load_info(info)
    # Read file info.dat
    file = string(info.folder,"/info.",info.ensemble,".dat")
    if isfile(file)
        fID = open(file,"r");   lines = readlines(fID);     close(fID)

        # Extract information
        pos1 = findfirst(": ",lines[1])
        info.moltype = lines[1][pos1[end]+1:end]
        pos2 = findfirst(": ",lines[2])
        info.dt = parse(Float64,lines[2][pos2[end]+1:end])
        pos3 = findfirst(": ",lines[3])
        info.natoms = parse(Int16,lines[3][pos3[end]+1:end])
    else error("File \"info.dat\" is empty") end

    return info
end

# Loading Thermo File
function load_thermo(info)
    # Read file info.dat and check for right format
    file = string(info.folder,"/thermo.",info.ensemble,".dat")
    if isfile(file)
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

        thermodat = thermo_dat(step, time, T, p, ρ, Etot, Ekin, Epot)
    else error("File \"thermo.dat\" is empty") end
end

# Loading Pressure File
function load_pressure(info)
    # Read file pressure.dat and check for right format
    file = string(info.folder,"/pressure.",info.ensemble,".dat")
    if isfile(file)
        fID = open(file,"r"); readline(fID); line2 = readline(fID); close(fID)
        if line2 != "# TimeStep v_V v_T c_thermo_press c_thermo_press[1] c_thermo_press[2] c_thermo_press[3] c_thermo_press[4] c_thermo_press[5] c_thermo_press[6]"
            error("Format of File \"pressure.dat\" not right")
        end

        # Read data
        dat = readdlm(file, skipstart=2)
        step = dat[:,1];    time = step*info.dt
        V    = dat[:,2];    T    = dat[:,3];    p    = dat[:,4]
        pxx  = dat[:,5];    pyy  = dat[:,6];    pzz  = dat[:,7]
        pxy  = dat[:,8];    pxz  = dat[:,9];    pyz  = dat[:,10]

        pdat = pressure_dat(step, time, V, T, p, pxx, pyy, pzz, pxy, pxz, pyz)
    else println("File \"pressure.dat\" is empty"); pdat = [] end
end

# Loading Dump File
function load_dump(info)
    # Initialization of posdat
    posdat = []

    # Read dump file timestep by timestep
    file = string(info.folder,"/atoms_position.",info.ensemble,".dump")
    if isfile(file)
        fID = open(file,"r")
        while !eof(fID)
            # Reading timestep and calculating time
            if (readline(fID) == "ITEM: TIMESTEP")
                step = parse(Int64,readline(fID))
                t = step*info.dt
            else error("File format wrong") end
            # Reading number of atoms
            if (readline(fID) == "ITEM: NUMBER OF ATOMS")
                natoms = parse(Int16,readline(fID))
            else error("File format wrong") end
            # Reading boundary coords
            if (readline(fID) == "ITEM: BOX BOUNDS pp pp pp")
                bounds = zeros(3,2)
                bounds[1,:] = parse.(Float64,split(readline(fID)))
                bounds[2,:] = parse.(Float64,split(readline(fID)))
                bounds[3,:] = parse.(Float64,split(readline(fID)))
            else error("File format wrong") end
            # Reading positions of atoms
            if (readline(fID) == "ITEM: ATOMS id mol xu yu zu ")
                id = zeros(natoms); molid = zeros(natoms);
                x = zeros(natoms); y = zeros(natoms); z = zeros(natoms);
                for i = 1:natoms
                    line_float = parse.(Float64,split(readline(fID)))
                    id[i] = line_float[1]; molid[i] = line_float[2]
                    x[i] = line_float[3]; y[i] = line_float[4]; z[i] = line_float[5]
                end
            else error("File format wrong") end

            posdat = [posdat;dump_dat(step, t, bounds, id, molid, x, y, z)]
        end
    else println("File \"atoms_positions.dump\" is empty") end
    return posdat
end

# Loading Heat Flux File
function load_heatflux(info)
    # Read file heat_flux.dat and check for right format
    file = string(info.folder,"/heat_flux.",info.ensemble,".dat")
    if isfile(file)
        fID = open(file,"r"); readline(fID); line2 = readline(fID); close(fID)
        if line2 != "# TimeStep v_V v_T c_j[1] c_j[2] c_j[3]"
            error("Format of File \"heat_flux.dat\" not right")
        end

        # Read data
        dat = readdlm(file, skipstart=2)
        step = dat[:,1];    time = step*info.dt
        V    = dat[:,2];    T    = dat[:,3]
        jx   = dat[:,4];    jy   = dat[:,5];    jz = dat[:,6]

        jdat = heat_dat(step, time, V, T, jx, jy, jz)
    else
        println("File \"heat_flux.dat\" is empty"); jdat = []
    end
end

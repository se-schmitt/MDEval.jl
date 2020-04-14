# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - Init
# Initializaiton of used structures and general functions
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

# Load Packages
using Dates
using DelimitedFiles
using Distributed
using Plots
using Printf
using Statistics

if no_procs > 1
    addprocs(no_procs)
elseif no_procs == 1 && nprocs()>1
    rmprocs(workers())
end

@everywhere using FFTW
@everywhere using StatsBase

# Physical Constants
global kB = 1.380649e-23

# Structures
# Structure to store general evaluation and simulation options
mutable struct info_struct
    folder
    ensemble
    do_multi
    n_equ
    moltype
    dt
    natoms
end

# Data structure to store thermo data (thermo.dat)
mutable struct thermo_dat
    step
    t
    T
    p
    ρ
    Etot
    Ekin
    Epot
end

# Data structure to store pressure tensor (pressure.dat)
mutable struct pressure_dat
    step
    t
    V
    T
    p
    pxx
    pyy
    pzz
    pxy
    pxz
    pyz
end

# Data structure to store atoms positions
mutable struct dump_dat
    step::Int64
    t::Float64
    bounds
    id
    molid
    x
    y
    z
end

# Data structure to store heat flux vector (heat_flux.dat)
mutable struct heat_dat
    step
    t
    V
    T
    jx
    jy
    jz
end

# Data structure to store single data
mutable struct single_dat
    val
    std
    err
end

# Data structure to store results
mutable struct results_struct
    T
    p
    ρ
    Etot
    Ekin
    Epot
    η
    η_V
    D
    λ
end

# Functions
# Function to convert to or froam DATA
function DATA2dats(DATA)
    thermodat   = DATA[1];  pdat    = DATA[2]
    posdat  = DATA[3];      jdat    = DATA[4]
    return thermodat, pdat, posdat, jdat
end
function dats2DATA(thermodat, pdat, posdat, jdat)
    DATA = Array{Any,1}(undef,4)
    DATA[1] = thermodat;    DATA[2] = pdat
    DATA[3] = posdat;       DATA[4] = jdat
    return DATA
end

# Function to get all subfolder
function get_subfolder(folder)
    # Get all subfolders
    paths = string.(folder,"/",readdir(folder))
    what = isdir.(paths)

    subfolder = paths[what]
end

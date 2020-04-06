# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - Init
# Initializaiton of used structures and general functions
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

# Load Packages
using Dates
using DelimitedFiles
using Plots
using Printf
using Statistics
using StatsBase

# Physical Constants
global kB = 1.380649e-23

# Structures
# Structure to store general evaluation and simulation options
mutable struct info_struct
    folder
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

# Function to integrate by trapezodial rule
function trapz(x,y,cum)
    # cum: 0 - output just total integral as float,
    #      1 - output vector with cumulative calculated integal
    # Check input
    if (length(x) != size(y,1)) error("Length of input vectors not equal") end

    # Numerical integration by trapeziodal rule
    if cum != 1
        dx = x[2:end] - x[1:end-1]
        int = sum(dx.*y[1:end-1,:] .+ (y[2:end,:]-y[1:end-1,:]).*dx./2, dims=1)
    else
        int = 0;
        for i=2:length(x) int = vcat(int,trapz(x[1:i],y[1:i,:],0)) end
    end
    return int
end

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
using Distributions
using Printf
using Statistics

if (nprocs() > 1) rmprocs(workers()) end
if (no_procs > 1) addprocs(no_procs) end

@everywhere using FFTW
@everywhere using LsqFit
@everywhere using Plots
@everywhere using StatsBase

# Physical Constants
global kB = 1.380649e-23    # J/K
global NA = 6.02214076e23   # 1/mol

# Structures
# Structure to store general evaluation and simulation options
mutable struct info_struct
    folder
    ensemble
    n_equ
    moltype
    dt
    natoms
    molmass
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
    id::Array{Int64,1}
    molid::Array{Int64,1}
    mass
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

# Data strucutre to store TDM settings
@everywhere mutable struct set_TDM
    folder::String
    subfolder::Array{String,1}
    do_out::Bool
    tskip::Float64
    cutcrit::Float64
    tcut
    name::String
    symbol::String
    unit::String
    nboot::Int64
end

# Data structure to save state information
mutable struct state_info
    T::Float64
    p::Float64
    ρ::Float64
    n::Float64
end

# Functions

# Function to get all subfolder
function get_subfolder(folder)
    # Get all subfolders
    paths = string.(folder,"/",readdir(folder))
    what = isdir.(paths) .& .!(occursin.("TransportProperties",paths))
    subf = paths[what]
end

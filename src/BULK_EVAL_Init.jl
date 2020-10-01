# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - Init
# Initializaiton of used structures and general functions
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

# Load Packages
using Dates
using DelimitedFiles
using Dierckx
using Distributed
using Distributions
using NLsolve
using Printf
using Statistics

if (nprocs() > 1) rmprocs(workers()) end
if (no_procs > 1) addprocs(no_procs) end

@everywhere using FFTW
@everywhere using LsqFit
# @everywhere using Plots
@everywhere using PyCall
@everywhere pygui(:qt)
@everywhere using PyPlot
close("all")
@everywhere using StatsBase

# Physical Constants
global kB = 1.380649e-23        # J/K
global NA = 6.02214076e23       # 1/mol
global eV2J = 1.602176565e-19   # 1 eV = 1.602176565e-19 J

# Structures
# Structure to store input variables read in from 'INPUT.txt'
mutable struct input_struct
    modus::String
    folders::Array{String,1}
    ensemble::String
    n_equ::Int64
    do_eval::Int64
    do_state::Int64
    n_boot::Int64
end

# Structure to store general evaluation and simulation options
mutable struct info_struct
    folder::String
    ensemble::String
    n_equ::Int64
    moltype::String
    dt::Float64
    natoms::Int64
    molmass::Float64
end

# Data structure to store thermo data (thermo.dat)
mutable struct thermo_dat
    step::Array{Int64,1}
    t::Array{Float64,1}
    T::Array{Float64,1}
    p::Array{Float64,1}
    ρ::Array{Float64,1}
    Etot::Array{Float64,1}
    Ekin::Array{Float64,1}
    Epot::Array{Float64,1}
end

# Data structure to store pressure tensor (pressure.dat)
mutable struct pressure_dat
    step::Array{Int64,1}
    t::Array{Float64,1}
    V::Array{Float64,1}
    T::Array{Float64,1}
    p::Array{Float64,1}
    pxx::Array{Float64,1}
    pyy::Array{Float64,1}
    pzz::Array{Float64,1}
    pxy::Array{Float64,1}
    pxz::Array{Float64,1}
    pyz::Array{Float64,1}
end

# Data structure to store atoms positions
mutable struct dump_dat
    step::Int64
    t::Float64
    bounds::Array{Float64,2}
    id::Array{Int64,1}
    molid::Array{Int64,1}
    type::Array{Int64,1}
    mass::Array{Float64,1}
    x::Array{Float64,1}
    y::Array{Float64,1}
    z::Array{Float64,1}
end

# Data structure to store heat flux vector (heat_flux.dat)
mutable struct heat_dat
    step::Array{Int64,1}
    t::Array{Float64,1}
    V::Array{Float64,1}
    T::Array{Float64,1}
    jx::Array{Float64,1}
    jy::Array{Float64,1}
    jz::Array{Float64,1}
end

# Data structure to store single data
mutable struct single_dat
    val::Float64
    std::Float64
    err::Float64
end

# Data structure to store results
mutable struct results_struct
    T
    p
    ρ
    Etot
    Ekin
    Epot
    c
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
    tcut::Float64
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
    n::Int64
    m::Float64
end

# Data structure to save profile information
mutable struct profile_data
    timestep::Array{Float64,1}
    id_chunk::Array{Float64,2}
    x::Array{Float64,2}
    Ncount::Array{Float64,2}
    ρn::Array{Float64,2}
    ρm::Array{Float64,2}
    T::Array{Float64,2}
    pxx::Array{Float64,2}
    pyy::Array{Float64,2}
    pzz::Array{Float64,2}
    pxy::Array{Float64,2}
    pxz::Array{Float64,2}
    pyz::Array{Float64,2}
end

# Functions
# Function to get all subfolders
function get_subfolder(folder)
    # Get all subfolders
    paths = string.(folder,"/",readdir(folder))
    what = isdir.(paths) .& .!(occursin.("TransportProperties",paths))
    subf = paths[what]
end

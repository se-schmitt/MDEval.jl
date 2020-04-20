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
global kB = 1.380649e-23

# Structures
# Structure to store general evaluation and simulation options
mutable struct info_struct
    folder
    ensemble
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
    what = isdir.(paths) .& .!(occursin.("TransportProperties",paths))
    subfolder = paths[what]
end

# Average of T, p, ρ
function ave_state(subf)
    # Make array
    n = length(subf)
    Tmat = zeros(n)
    pmat = zeros(n)
    ρmat = zeros(n)
    for i = 1:n
        res = load_result(string(subfolder[i],"/result.dat"))
        Tmat[i] = res.T.val
        pmat[i] = res.p.val
        ρmat[i] = res.ρ.val
    end

    # Save in structure
    T = single_dat(mean(Tmat), std(Tmat), std(Tmat)/n)
    p = single_dat(mean(pmat), std(pmat), std(pmat)/n)
    ρ = single_dat(mean(ρmat), std(ρmat), std(ρmat)/n)

    return T, p, ρ
end

# Load result file
function load_result(file)
    fID = open(file,"r")
    lines = readlines(fID);
    close(fID)
    res = results_struct([],[],[],[],[],[],[],[],[],[])

    for i = 1:length(lines)
        pos_colon = findfirst(isequal(':'),lines[i])
        pos_open = findfirst(isequal('('),lines[i])
        pos_comma = findfirst(isequal(','),lines[i])
        pos_close = findfirst(isequal(')'),lines[i])
        if length(lines[i])>1 && !(lines[i][1] == '#')
            if !occursin("---",lines[i])
                if lines[i][1] == 'ρ' || lines[i][1] == 'η'
                    name = lines[i][1:pos_colon-2]
                else
                    name = lines[i][1:pos_colon-1]
                end
                if isnothing(pos_open) && isnothing(pos_close)
                    val = parse(Float64,strip(lines[i][pos_colon+1:pos_comma-1]))
                    std = []; err = []
                elseif pos_close < pos_comma
                    val = parse(Float64,strip(lines[i][pos_colon+1:pos_open-1]))
                    std = parse(Float64,strip(lines[i][pos_open+1:pos_close-1]))
                    err = []
                elseif (pos_open < pos_comma) && (pos_close > pos_comma)
                    val = parse(Float64,strip(lines[i][pos_colon+1:pos_open-1]))
                    std = parse(Float64,strip(lines[i][pos_open+1:pos_comma-1]))
                    err = parse(Float64,strip(lines[i][pos_comma+1:pos_close-1]))
                end
                if (name == "T") res.T = single_dat(val,std,err) end
                if (name == "p") res.p = single_dat(val,std,err) end
                if (name == "ρ") res.ρ = single_dat(val,std,err) end
                if (name == "Etot") res.Etot = single_dat(val,std,err) end
                if (name == "Ekin") res.Ekin = single_dat(val,std,err) end
                if (name == "Epot") res.Epot = single_dat(val,std,err) end
                if (name == "η") res.η = single_dat(val,std,err) end
                if (name == "η_") res.η_V = single_dat(val,std,err) end
                if (name == "D") res.D = single_dat(val,std,err) end
                if (name == "λ") res.λ = single_dat(val,std,err) end
            end
        end
    end
    return res
end

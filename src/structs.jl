# Structure to store general evaluation and simulation options
mutable struct Info
    folder::String
    ensemble::String
    n_equ::Int64
    moltype::String
    dt::Float64
    natoms::Int64
    molmass::Float64
    n_blocks::Int64
    N_bin::Union{Int64,Missing}
    r_cut::Union{Float64,Missing}
end

# Data structure to store thermo data (thermo.dat)
mutable struct ThermoDat
    step::Array{Int64,1}
    t::Array{Float64,1}
    T::Array{Float64,1}
    p::Array{Float64,1}
    ρ::Array{Float64,1}
    Etot::Array{Float64,1}
    Ekin::Array{Float64,1}
    Epot::Array{Float64,1}
    pyz::Array{Float64,1}
    Qhot::Array{Float64,1}
    Qcold::Array{Float64,1}
    method::String
end

# Data structure to store pressure tensor (pressure.dat)
mutable struct PressureDat
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
mutable struct DumpDat
    step::Int64
    t::Float64
    bounds::Array{Float64,2}
    id::Array{Int64,1}
    type::Array{Int64,1}
    molid::Array{Int64,1}
    moltype::Array{Int64,1}
    mass::Array{Float64,1}
    x::Array{Float64,1}
    y::Array{Float64,1}
    z::Array{Float64,1}
end

# Data structure to store heat flux vector (heat_flux.dat)
mutable struct HeatDat
    step::Array{Int64,1}
    t::Array{Float64,1}
    V::Array{Float64,1}
    T::Array{Float64,1}
    jx::Array{Float64,1}
    jy::Array{Float64,1}
    jz::Array{Float64,1}
end

# Data structure to store single data
mutable struct SingleDat
    val::Float64
    std::Float64
    err::Float64
end

# Data structure to store results
mutable struct ResultsDat
    T
    p
    ρ
    x
    Etot
    Ekin
    Epot
    c
    η
    η_V
    D
    λ
    Rg
end

mutable struct ResultsDatNEMD
    T
    p
    ρ
    x
    Etot
    Ekin
    Epot
    η
    s_rate
    λ
end

# Data strucutre to store TDM settings
mutable struct OptsTDM
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

# Data structure to store thermo data (thermo.vle.2phase.dat)
mutable struct ThermoVLEDat
    step::Array{Int64,1}
    t::Array{Float64,1}
    T::Array{Float64,1}
    px::Array{Float64,1}
    py::Array{Float64,1}
    pz::Array{Float64,1}
    ρ::Array{Float64,1}
    Etot::Array{Float64,1}
    Ekin::Array{Float64,1}
    Epot::Array{Float64,1}
end

# Data structure to save profile information of vle simulations
mutable struct ProfileVLEDat
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

# Data structure to save profile information of shear simulations
mutable struct ProfileShearDat
    timestep::Array{Float64,1}
    id_chunk::Array{Float64,2}
    x::Array{Float64,2}
    Ncount::Array{Float64,2}
    vy::Array{Float64,2}
end

# Data structure to save profile information of heat simulations
mutable struct ProfileHeatDat
    timestep::Array{Float64,1}
    id_chunk::Array{Float64,2}
    x::Array{Float64,2}
    Ncount::Array{Float64,2}
    T::Array{Float64,2}
end

# Data structure to save profile information of rdf
mutable struct ProfileRDFDat
    timestep::Array{Float64,1}
    id_chunk::Array{Float64,2}
    r::Array{Float64,2}
    g_r::Array{Any,1}
    N_coord::Array{Any,1}
end


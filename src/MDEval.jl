module MDEval

# Imports
using Dates, Printf, ProgressMeter
using DelimitedFiles, Distributed
using FFTW, LsqFit, NLsolve
using Distributions, Statistics, StatsBase, RollingFunctions
using PyCall, PyPlot

# Globals
dateformat = "yyyy-mm-dd HH:MM:SS"
reduced_units = false
kB = 1.380649e-23                   # J/K
NA = 6.02214076e23                  # 1/mol
eV2J = 1.602176565e-19              # 1 eV = 1.602176565e-19 J

# Include files
include("structs.jl")
include("input.jl")
include("output.jl")
include("load.jl")
include("utils.jl")
include("emd_single.jl")
include("emd_state.jl")
include("nemd_shear.jl")
include("nemd_shear_state.jl")
include("nemd_heat.jl")
include("eval_vle.jl")
include("eval_structure.jl")
include("green_kubo.jl")
include("stokes_einstein.jl")
include("transport_properties.jl")

# Main function
"""
    mdeval(mode::Symbol, folder::String, keywords::NamedTuple)

Evaluate molecular dynamics (MD) simulation files.

# Arguments
- `mode::Symbol`: Mode of evaluation (`:single_run`, `:tdm`, `:vle`, `:nemd_shear`, `:nemd_heat`).
- `folder::String`: Folder containing the MD simulation files.
- `keywords::NamedTuple`: Keywords for evaluation.
"""
function mdeval(mode::Symbol, folder::String, keywords::NamedTuple)
    println("="^60,"\n","START: \t",Dates.format(now(),dateformat),"\nFolder: ",folder,"\n","-"^60)
    println("$(Dates.format(now(),dateformat)): RUNNING ... ")

    # Get evaluation options
    opts = get_opts(mode,folder,keywords)

    # Run evaluation
    if mode == :single_run
        eval_single(folder,opts)

    elseif mode == :tdm
        # Get all subfolder
        subfolder = get_subfolder(folder)
        if isempty(subfolder)
            subfolder = [folder]
        end

        # Evaluate all subfolders
        if opts.do_single
            eval_single.(subfolder,Ref(opts))
        end

        # Evaluate state folder
        if opts.do_state
            eval_state(subfolder, opts)
        end

    elseif mode == :vle
        eval_vle(folder,opts)

    elseif mode == :nemd_shear
        # Get all subfolder
        subfolder = get_subfolder(folder)
        if isempty(subfolder)
            subfolder = [folder]
        end

        # Evaluate all subfolders
        if opts.do_single
            eval_shear.(subfolder,Ref(opts))
        end

        # Evaluate state folder
        if opts.do_state && length(subfolder) > 1
            eval_state_shear(subfolder, opts)
        end

    elseif mode == :nemd_heat
        eval_heat(folder,opts)
    end
    close("all")
    
    println("$(Dates.format(now(),dateformat)): DONE")
    println("="^60)
end 

# Export 
export mdeval

end

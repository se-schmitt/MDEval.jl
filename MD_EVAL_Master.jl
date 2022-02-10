## MD_EVAL_Master.jl
# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - Master
# Calculation of Averaged Thermodynamic Properties, Transport Coefficients, etc.
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

# --- Global Settings ---
# Number of processors to use
if length(ARGS) >= 2
    no_procs = parse(Int64,ARGS[2])
else
    no_procs = 1
end

# --- Define Functions and Structures ---
include("src/Init.jl")
include("src/Main.jl")
include("src/LoadData.jl")
include("src/EvalSingle.jl")
include("src/EvalState.jl")
include("src/EvalVLE.jl")
include("src/OutputResult.jl")
include("src/GreenKubo.jl")
include("src/StokesEinstein.jl")
include("src/TransportProperties.jl")
include("src/EvalStructure.jl")
include("src/EvalNEMDShear.jl")
include("src/EvalStateNEMDShear.jl")
include("src/EvalNEMDHeat.jl")
# --- Start Program ---
main(ARGS)

## MD_EVAL_Master.jl
# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - Master
# Calculation of Averaged Thermodynamic Properties, Transport Coefficients, etc.
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

# --- Global Settings ---
# Number of processors to use
no_procs = 1

# --- Define Functions and Structures ---
include("src/MD_EVAL_Init.jl")
include("src/MD_EVAL_Main.jl")
include("src/MD_EVAL_LoadData.jl")
include("src/MD_EVAL_EvalSingle.jl")
include("src/MD_EVAL_EvalState.jl")
include("src/MD_EVAL_EvalVLE.jl")
include("src/MD_EVAL_OutputResult.jl")
include("src/MD_EVAL_GreenKubo.jl")
include("src/MD_EVAL_StokesEinstein.jl")
include("src/MD_EVAL_TransportProperties.jl")
include("src/MD_EVAL_EvalStructure.jl")

# --- Start Program ---
main()

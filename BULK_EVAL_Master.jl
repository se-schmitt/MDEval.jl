# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - Master
# Calculation of Averaged Thermodynamic Properties, Transport Coefficients, etc.
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

# --- Global Settings ---
# Number of processors to use
no_procs = 4
# Calculating in reduced units
reduced_units = true

# --- Define Functions and Structures ---
include("src/BULK_EVAL_Init.jl")
include("src/BULK_EVAL_Main.jl")
include("src/BULK_EVAL_LoadData.jl")
include("src/BULK_EVAL_EvalSingle.jl")
include("src/BULK_EVAL_EvalState.jl")
include("src/BULK_EVAL_OutputResult.jl")
include("src/BULK_EVAL_TransportProperties.jl")

# --- Start Program ---
main()

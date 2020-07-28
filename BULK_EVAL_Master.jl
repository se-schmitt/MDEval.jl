# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - Master
# Calculation of Averaged Thermodynamic Properties, Transport Coefficients, etc.
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

no_procs = 2

reduced_units = true

include("src/BULK_EVAL_Init.jl")
include("src/BULK_EVAL_Main.jl")
include("src/BULK_EVAL_LoadData.jl")
include("src/BULK_EVAL_EvalData.jl")
include("src/BULK_EVAL_EvalState.jl")
include("src/BULK_EVAL_OutputResult.jl")
include("src/BULK_EVAL_TransportProperties.jl")

main()

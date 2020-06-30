# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - Master
# Calculation of Averaged Thermodynamic Properties, Transport Coefficients, etc.
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

no_procs = 1

include("Software/BULK_EVAL_Init.jl")
include("Software/BULK_EVAL_Main.jl")
include("Software/BULK_EVAL_LoadData.jl")
include("Software/BULK_EVAL_EvalData.jl")
include("Software/BULK_EVAL_EvalState.jl")
include("Software/BULK_EVAL_OutputResult.jl")
include("Software/BULK_EVAL_TransportProperties.jl")

main()

# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - Master
# Calculation of Averaged Thermodynamic Properties, Transport Coefficients, etc.
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

no_procs = 8            # Parallel computing if no_cpu_procs > 1

include("Software/BULK_EVAL_Init.jl")
include("Software/BULK_EVAL_LoadData.jl")
include("Software/BULK_EVAL_EvalData.jl")
include("Software/BULK_EVAL_OutputResult.jl")
include("Software/BULK_EVAL_TransportProperties.jl")
println(string("START: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS"),"\n---"))

# ---------- INPUT ----------
# Data to read
folder = "F:/MD_Bulk/Methane/Test_multiFiles"
ensemble = "NVT"
# Evaluation
n_equ = 0                   # Timesteps to wait from start of simulations
do_eval = 1                 # set 0, if single folders already evaluated
do_state = 0                # set 1, if get averaged values of subfolder
nboot = 200                 # Number of bootstrapping repetitions
# ---------------------------

# Get all subfolder
subfolder = get_subfolder(folder)
if isempty(subfolder)
    subfolder = [folder]
end

# Evaluation of single folder
if do_eval == 1
    for i = 1:length(subfolder)
        ifolder = subfolder[i]
        # Initialization of info structure
        info = info_struct(ifolder,ensemble,n_equ,"",[],[])

        # Load data
        data, info = LoadData(info)
        println(string("LoadData DONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS")))

        # Evaluate Data
        results = EvalData(data, info)
        println(string("EvalData DONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS")))

        # Output Results
        OutputResult(results, ifolder)

        println(string("---\nFolder ",i," / ",length(subfolder)," DONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS"),"\n---"))
    end
end

# Averaging simulations
if do_state == 1
    # T, p, ρ
    T, p, ρ = ave_state(subfolder)
    println(string("ave_state DONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS")))

    # Transport Properties
    η, λ = TransportProperties(set_TDM(folder,subfolder,false,2.0,0.4,[],"","","",nboot))
    println(string("TransportProperties DONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS")))
end

# Output Data

# END
if nprocs() > 1
    rmprocs(workers())
end
println(string("---\nDONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS")))

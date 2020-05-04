# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - Master
# Calculation of Averaged Thermodynamic Properties, Transport Coefficients, etc.
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

# ---------- INPUT ---------
# Data to read
folders = ["F:/MD_Bulk/Decane/Study_ForceFields/SIM_T_298.15K_p_0.1MPa"]
ensemble = "NVT"
# Evaluation
n_equ = 0                   # Timesteps to wait from start of simulations
do_eval = 1                 # set 0, if single folders already evaluated
do_state = 1                # set 1, if get averaged values of subfolder
nboot = 5                   # Number of bootstrapping repetitions
# Settings
no_procs = 4                # Parallel computing if no_cpu_procs > 1
# ---------------------------

include("Software/BULK_EVAL_Init.jl")
include("Software/BULK_EVAL_LoadData.jl")
include("Software/BULK_EVAL_EvalData.jl")
include("Software/BULK_EVAL_EvalState.jl")
include("Software/BULK_EVAL_OutputResult.jl")
include("Software/BULK_EVAL_TransportProperties.jl")
println(string("START: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS"),"\n__________"))

for folder in folders
    # Get all subfolder
    subfolder = get_subfolder(folder)
    if isempty(subfolder)
        subfolder = [folder]
    end
    println("Folder: ",folder)

    # Evaluation of single folder
    if do_eval == 1
        for i = 1:length(subfolder)
            # Initialization of info structure
            info = info_struct(subfolder[i],ensemble,n_equ,"",NaN,0,NaN)

            # Evaluate Data
            EvalData(info)
            println(string("\n---Folder ",i," / ",length(subfolder)," DONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS")))
        end
        println("__________")
    end

    # Averaging simulations
    if do_state == 1
        EvalState(subfolder)
        println(string("\n---State Evaluation DONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS")))
        println("__________")
    end
end

# END
if nprocs() > 1
    rmprocs(workers())
end
println(string("---\nDONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS")))

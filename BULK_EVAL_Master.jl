# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - Master
# Calculation of Averaged Thermodynamic Properties, Transport Coefficients, etc.
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

no_procs = 6            # Parallel computing if no_cpu_procs > 1

include("Software/BULK_EVAL_Init.jl")
include("Software/BULK_EVAL_LoadData.jl")
include("Software/BULK_EVAL_EvalData.jl")
include("Software/BULK_EVAL_OutputResult.jl")
println(string("START: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS"),"\n---"))

# ---------- INPUT ----------
# Data to read
folder = "C:/Users/Schmitt/Desktop/TEST_EvalBulk/Test_compass"
ensemble = "NVT"
do_multi = 0                # 1 - Evaluate subfolder of folder
# Evaluation
n_equ = 0                   # Timesteps to wait from start of simulations
single_state = 1            # 1 - Simulations to Evaluate at Same State Point
do_TDM = 0                  # 1 - Evaluate by Time Decomposition Method ()
# ---------------------------

# Get all subfolder
if do_multi == 1
     subfolder = get_subfolder(folder); nfolder = size(subfolder)
else subfolder = [folder];                nfolder = 1 end

# Initialization of Results Array
RESULTS = Vector{Any}(undef,nfolder)

k = 0
for ifolder in subfolder
    global k += 1;

    # Initialization of info structure
    info = info_struct(ifolder,ensemble,do_multi,n_equ,"",[],[])

    # Load data
    DATA, info = LoadData(info)
    println(string("ave_thermo DONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS")))

    # Evaluate Data
    RESULTS[k] = EvalData(DATA, info)
    println(string("ave_thermo DONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS")))
    # alle zeitabhängigen Daten (wie η(t), acf(t), msd(t), D(t), λ(t), acf_λ(t))
    # werden als Textdateien ausschreiben (+ evtl. als Plot)

    # Output Results
    OutputResult(RESULTS[k], ifolder)

    println(string("---\nFolder ",k," / ",length(subfolder)," DONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS"),"\n---"))
end

# Time Decomposition Method
if do_TDM == 1
    TimeDecompositionMethod()
end

# Output Data

# END
if nprocs() > 1
    rmprocs(workers())
end
println(string("---\nDONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS")))

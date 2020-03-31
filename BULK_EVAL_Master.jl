# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - Master
# Calculation of Averaged Thermodynamic Properties, Transport Coefficients, etc.
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

include("Software/BULK_EVAL_Init.jl")
include("Software/BULK_EVAL_LoadData.jl")
include("Software/BULK_EVAL_EvalData.jl")

# ---------- INPUT ----------
# Data to read
folder = "C:/Daten/1_Projekte/MD_Bulk/Input_Skripte_EMD/MinBSP_Pentane"
do_multi = 0;           # 1 - Evaluate subfolder of folder
# Evaluation
# n_eval = 2e6;           # Timesteps to evaluate (length of evaluated sections)
# n_wait = 1e6;           # Timesteps to wait to create independent trajectories
n_equ = 1000;              # Timesteps to wait from start of simulations
# do_TDM = 0;             # 1 - Evaluate by Time Decomposition Method ()
# ---------------------------

# Get all subfolder
if do_multi == 1
     subfolder = get_subfolder(folder); nfolder = size(subfolder)
else subfolder = [folder];                nfolder = 1 end

# Initialization of Results Array
RESULTS = Vector{Any}(undef,nfolder)

for ifolder in subfolder
    # Initialization of info structure
    info = info_struct(ifolder,do_multi,n_equ,"",[],[])

    # Load data
    DATA, info = LoadData(info)

    # Evaluate Data
    RESULTS = EvalData(DATA, info)
    # alle Berechnungsdaten wie η(t), acf(t), msd(t), D(t), λ(t), acf_λ(t), ...
    # als Textdateien ausschreiben (+ evtl. als Plot) -> in RESULTS nur End-
    # ergebnisse speichern

end

# END
println("---")
println(string("DONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS")))

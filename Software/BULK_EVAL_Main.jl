# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - Main
# Contains main function to run evaluation of bulk simulations
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

# include("Software/BULK_EVAL_Init.jl")


# include("Software/BULK_EVAL_Init.jl")
# include("Software/BULK_EVAL_LoadData.jl")
# include("Software/BULK_EVAL_EvalData.jl")
# include("Software/BULK_EVAL_EvalState.jl")
# include("Software/BULK_EVAL_OutputResult.jl")
# include("Software/BULK_EVAL_TransportProperties.jl")

function main()
    println(string("----------\nSTART: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS"),"\n----------"))
    println(string("Number of processors: ",no_procs))

    # Read input parameters
    inpar = read_input()
    error()

    # Loop over all folders
    for folder in inpar.folders
        # Get all subfolder
        subfolder = get_subfolder(folder)
        if isempty(subfolder)
            subfolder = [folder]
        end
        println("Folder: ",folder)

        # Evaluation of single folder
        if inpar.do_eval == 1
            for i = 1:length(subfolder)
                # Initialization of info structure
                info = info_struct(subfolder[i],inpar.ensemble,inpar.n_equ,"",NaN,0,NaN)

                # Evaluate Data
                EvalData(info)
                println(string("\n---\nFolder ",i," / ",length(subfolder)," DONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS")))
            end
            println("----------")
        end

        # Averaging simulations
        if inpar.do_state == 1
            EvalState(subfolder)
            println(string("\n---State Evaluation DONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS")))
            println("----------")
        end
    end

    # END
    if nprocs() > 1
        rmprocs(workers())
    end
    println(string("----------\nDONE: ",Dates.format(now(),"yyyy-mm-dd HH:MM:SS"),"\n----------"))
end

# Subfunctions
# Function to read input file
function read_input()
    # Filename
    file = "INPUT.txt"

    # Initialization of input variables
    inpar = input_struct([],"",-1,-1,-1,-1)

    fID = open(file,"r")
    while !eof(fID)
        line = readline(fID)
        if line == "#folder"
            inpar.folders = vcat(inpar.folders,readline(fID))
        elseif line == "#ensemble"
            inpar.ensemble = readline(fID)
        elseif line == "#timesteps_EQU"
            inpar.n_equ = parse(Int64,readline(fID))
        elseif line == "#DO_evaluation"
            inpar.do_eval = parse(Int64,readline(fID))
        elseif line == "#DO_state"
            inpar.do_state = parse(Int64,readline(fID))
        elseif line == "#N_boot"
            inpar.n_boot = parse(Int64,readline(fID))
        end
    end

    # Check input structure
    if (inpar.n_equ == -1 | inpar.do_eval == -1 |
        inpar.do_state == -1 | inpar.n_boot == -1)
        error("-1 in Input Structure: Something is wrong!")
    end

    return inpar
end

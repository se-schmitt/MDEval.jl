# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - Main
# Contains main function to run evaluation of bulk simulations
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

function main()
    sline = "------------------------------------------------------------------"
    dline = "=================================================================="
    fdate = "yyyy-mm-dd HH:MM:SS"
    println(dline)
    println(string("START: \t",Dates.format(now(),fdate)))
    println(sline)
    println(string("Number of processors: ",no_procs))
    println(sline)

    # Read input parameters
    inpar = read_input()

    # Loop over all folders
    for folder in inpar.folders
        println(string("Folder: ",folder))

        ## Mode "single_run" ---------------------------------------------------
        if inpar.mode == "single_run"
            print(string(Dates.format(now(),fdate),": RUNNING ... "))

            # Evaluate Data
            EvalSingle(folder,inpar)

            println(string("   →    ",Dates.format(now(),fdate),": DONE"))

        ## Mode "tdm" ----------------------------------------------------------
        elseif inpar.mode == "tdm"
            # Get all subfolder
            subfolder = get_subfolder(folder)

            if isempty(subfolder)
                subfolder = [folder]
            end

            # Evaluation of single folder
            if inpar.do_eval == 1
                for i = 1:length(subfolder)
                    print(string("Subfolder ",i," / ",length(subfolder)," → RUNNING ... "))

                    # Evaluate Data
                    EvalSingle(subfolder[i],inpar)

                    println(string("   →    ",Dates.format(now(),fdate),": DONE"))
                end
                println(sline)
            end

            # Averaging simulations
            if inpar.do_state == 1
                println("State Evaluation")
                intro = string(Dates.format(now(),fdate),": RUNNING ... ")
                println(intro)

                # Do state evaluation
                EvalState(subfolder, inpar)

                println(string(intro,"   →    ",Dates.format(now(),fdate),": DONE"))
            end

        ## Mode "vle" ----------------------------------------------------------
        elseif inpar.mode == "vle"
            info = info_struct(folder,inpar.ensemble,inpar.n_equ,"",NaN,0,NaN)
            EvalVLE(info)
        end
        println(sline)
    end

    # Output in single file
    if inpar.mode == "transport"
        dlm_output(inpar.folders)
    elseif inpar.mode == "vle"
    end

    # END
    if nprocs() > 1
        rmprocs(workers())
    end

    println(string(Dates.format(now(),fdate),": DONE"))
    println(dline)
end

## Subfunctions
# Function to read input file
function read_input()
    # Filename
    file = "INPUT.txt"

    # Initialization of input variables
    inpar = input_struct("",[],"",-1,-1,-1,-1,-1,-1)

    fID = open(file,"r")
    while !eof(fID)
        line = readline(fID)
        if line == "#mode"
            inpar.mode = readline(fID)
        elseif line == "#folder"
            foldername = readline(fID)
            if foldername[end] == '*'
                pos = findlast(isequal('/'),foldername)
                folders = string.(foldername[1:pos],readdir(foldername[1:pos]))
                folders = folders[isdir.(folders)]

                if pos != length(foldername)-1
                    startstr = foldername[pos+1:end-1]
                    what = ones(length(folders)) .== zeros(length(folders))
                    for i = 1:length(folders)
                        f = folders[i]
                        what[i] = length(f) > pos+length(startstr) && f[pos+1:pos+length(startstr)] == startstr
                    end
                    folders = folders[what]
                end
                inpar.folders = vcat(inpar.folders,folders)
            else
                inpar.folders = vcat(inpar.folders,foldername)
            end
        elseif line == "#ensemble"
            inpar.ensemble = readline(fID)
        elseif line == "#timesteps_EQU"
            inpar.n_equ = parse(Int64,readline(fID))
        elseif line == "#DO_single"
            inpar.do_eval = parse(Int64,readline(fID))
        elseif line == "#DO_state"
            inpar.do_state = parse(Int64,readline(fID))
        elseif line == "#N_boot"
            inpar.n_boot = parse(Int64,readline(fID))
        elseif line == "#corr_length"
            inpar.corr_length = parse(Int64,readline(fID))
        elseif line == "#span_corr_fun"
            inpar.span_corr_fun = parse(Int64,readline(fID))
        end
    end

    # Check input structure
    # Mode "single_run"
    if (inpar.n_equ == -1 || inpar.corr_length == -1 ||
        inpar.span_corr_fun == -1) && inpar.mode == "single_run"
        error("Parameters missing for mode \"single_run\"!")
    end
    # Mode "tdm"
    if (inpar.n_equ == -1 || inpar.do_eval == -1 ||
        inpar.do_state == -1 || inpar.n_boot == -1) && inpar.mode == "tdm"
        error("Parameters missing for mode \"tdm\"!")
    end

    return inpar
end

# Output results in dlm file
function dlm_output(folders)
    # Settings
    error_type = :std                       # :std, :err or :none
    props = [:T,:p,:ρ,:x,:η,:D,:λ]          # Properties to write out

    # Get results
    results = []
    folder_str = "Simulation folder:\n"
    for f in folders
        if isfile(string(f,"/result.dat"))
            results = vcat(results,load_result(string(f,"/result.dat")))
            folder_str = string(folder_str,f,"\n")
        else
            println(string("No result file in folder: \"",f,"\""))
        end
    end

    if !isempty(results)
        # Create matrix for values with error (stds/errs)
        Nmol = maximum(length.(getfield.(results,:x)))
        Nval = length(props) + 2 * (Nmol - 1)
        if error_type != :none
            mat = Array{Float64,2}(undef,length(results),Nval*2).*NaN
        else
            mat = Array{Float64,2}(undef,length(results),Nval).*NaN
        end

        # Prepare header
        vheader = []
        for i = 1:length(props)
            if props[i] in [:x,:D]
                for j = 1:Nmol
                    vheader = vcat(vheader,string(props[i],j))
                    if (error_type != :none) vheader = vcat(vheader,"+-") end
                end
            else
                vheader = vcat(vheader,string(props[i]))
                if (error_type != :none) vheader = vcat(vheader,"+-") end
            end
        end

        for i = 1:length(results)
            jval = 0
            for j = 1:length(props)
                sdat = getfield(results[i],props[j])

                if typeof(sdat) == single_dat
                    jval += 1
                    if error_type != :none
                        if typeof(sdat) == single_dat
                            mat[i,(jval-1)*2+1:jval*2] = [sdat.val,getfield(sdat,error_type)]
                        else
                            mat[i,(jval-1)*2+1:jval*2] = [NaN,NaN]
                        end
                    else
                        if typeof(sdat) == single_dat
                            mat[i,jval] = sdat.val
                        else
                            mat[i,jval] = NaN
                        end
                    end
                elseif typeof(sdat) in [Array{Any,1},Array{single_dat,1}] && length(sdat) > 0
                    for k = 1:Nmol
                        jval += 1
                        if k <= length(sdat)
                            if error_type != :none
                                if typeof(sdat[k]) == single_dat
                                    mat[i,(jval-1)*2+1:jval*2] = [sdat[k].val,getfield(sdat[k],error_type)]
                                else
                                    mat[i,(jval-1)*2+1:jval*2] = [NaN,NaN]
                                end
                            else
                                if typeof(sdat[k]) == single_dat
                                    mat[i,jval] = sdat[k].val
                                else
                                    mat[i,jval] = NaN
                                end
                            end
                        end
                    end
                else
                    jval += 1
                end
            end
        end

        # Delete columns with NaN
        nonan_col = (sum(isnan.(mat),dims=1) .!= size(mat,1))[:]
        mat = mat[:,nonan_col]
        vheader = vheader[nonan_col]

        header = ""
        for hstr in vheader
            header = string(header,hstr," ")
        end

        # Write Files
        f1 = folders[1]
        posf = findall(isequal('/'),f1)
        folder_save = ""
        for posTMP in posf
            inmat = findfirst.(f1[1:posTMP],folders)
            if sum(isnothing.(inmat)) == 0
                folder_save = f1[1:posTMP]
            else
                break
            end
        end
        datestr = Dates.format(now(),"yyyy-mm-dd_HH.MM")
        file_save = string(folder_save,"results_",datestr,".dat")
        println("Results saved in file: ",file_save)

        fID = open(file_save,"w")
        println(fID,header)
        # writedlm(fID,mat,' ')
        for i = 1:size(mat,1)
            for j = 1:size(mat,2)
                if vheader[j][1] == 'x'
                    @printf(fID,"%.3f ",mat[i,j])
                else
                    @printf(fID,"%.8e ",mat[i,j])
                end
            end
            @printf(fID,"\n")
        end
        close(fID)

        # Save file with folder names
        fID = open(string(folder_save,"results_",datestr,"_folder.dat"),"w")
        println(fID,folder_str)
        close(fID)
    end
end

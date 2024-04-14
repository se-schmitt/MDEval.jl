# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - EvalData
# Function to evaluate data from single MD run
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

# Function to ouput results of a single simulation / state point
function output_results(result::ResultsDat, folder::String)
    # Create Filepath
    path = string(folder,"/result.dat")

    no = 0
    while isfile(path)
        no += 1
        if no == 1  path = replace(path,".dat" => "_$(no).dat")
        else        path = replace(path,"_$(no-1).dat" => "_$(no).dat")  
        end
    end

    # Write to file
    fID = open(path,"w")
    line1 = "# Created by MD - Bulk Evaluation, folder: $folder, time: $(Dates.format(now(),"yyyy-mm-dd HH:MM:SS"))"
    line2 = "# Format: val (std dev., std err.); [T]=K, [p]=MPa, [ρ]=g/ml, [E]=eV, [η]=Pa*s, [D]=m^2/s, [λ]=W/(m*K)"
    if (reduced_units) line2 = "# Format: val (std dev., std err.); reduced units" end
    header = string(line1,"\n",line2,"\n")
    print(fID,header)

    # Write data
    print_prop(fID, result.T, "T")
    print_prop(fID, result.p, "p")
    print_prop(fID, result.ρ, "ρ")
    print_prop(fID, result.x, "x")
    print_prop(fID, result.Etot, "Etot")
    print_prop(fID, result.Ekin, "Ekin")
    print_prop(fID, result.Epot, "Epot")
    print_prop(fID, result.c, "c")
    print_prop(fID, result.η, "η")
    print_prop(fID, result.η_V, "η_V")
    print_prop(fID, result.D, "D")
    print_prop(fID, result.λ, "λ")
    print_prop(fID, result.Rg, "Rg")
    close(fID)
end

function output_resultsNEMD(result::ResultsDatNEMD, folder::String; note="")
    # Create Filepath
    path = string(folder,"/result.dat")

    # Write to file
    fID = open(path,"w")
    line1 = "# Created by MD - Bulk Evaluation, folder: $folder, time: $(Dates.format(now(),"yyyy-mm-dd HH:MM:SS"))"
    line2 = "# Format: val (std dev., std err.); [T]=K, [p]=MPa, [ρ]=g/ml, [E]=eV, [η]=Pa*s, [s_rate]=1/s"
    if (reduced_units) line2 = "# Format: val (std dev., std err.); reduced units" end
    header = string(line1,"\n",line2,note,"\n")
    print(fID,header)
    # Write data
    print_prop(fID, result.T, "T")
    print_prop(fID, result.p, "p")
    print_prop(fID, result.ρ, "ρ")
    print_prop(fID, result.x, "x")
    print_prop(fID, result.Etot, "Etot")
    print_prop(fID, result.Ekin, "Ekin")
    print_prop(fID, result.Epot, "Epot")
    print_prop(fID, result.η, "η")
    print_prop(fID, result.s_rate, "s_rate")
    print_prop(fID, result.λ, "λ")
    close(fID)
end

# Function to ouput data from vle simulations (single simulation)
function output_results_VLE(dat::ThermoVLEDat, folder::String)
    # Create Filepath
    path = string(folder,"/result.dat")

    # Write to file
    fID = open(path,"w")
    line1 = "# Created by MD - Bulk Evaluation, folder: $folder, time: $(Dates.format(now(),"yyyy-mm-dd HH:MM:SS"))"
    line2 = "# Format: val (std dev., std err.); [T]=K, [p]=MPa, [ρ]=g/ml, [E]=eV"
    if (reduced_units) line2 = "# Format: val (std dev., std err.); reduced units" end
    header = string(line1,"\n",line2,"\n")
    print(fID,header)

    # Calculation of means, stds, errs
    T = SingleDat(mean(dat.T),block_average(dat.T)[1],block_average(dat.T)[2])
    px = SingleDat(mean(dat.px),block_average(dat.px)[1],block_average(dat.px)[2])
    py = SingleDat(mean(dat.py),block_average(dat.py)[1],block_average(dat.py)[2])
    pz = SingleDat(mean(dat.pz),block_average(dat.pz)[1],block_average(dat.pz)[2])
    ρ = SingleDat(mean(dat.ρ),block_average(dat.ρ)[1],block_average(dat.ρ)[2])
    Etot = SingleDat(mean(dat.Etot),block_average(dat.Etot)[1],block_average(dat.Etot)[2])
    Ekin = SingleDat(mean(dat.Ekin),block_average(dat.Ekin)[1],block_average(dat.Ekin)[2])
    Epot = SingleDat(mean(dat.Epot),block_average(dat.Epot)[1],block_average(dat.Epot)[2])

    # Write data
    print_prop(fID, T, "T")
    print_prop(fID, px, "px")
    print_prop(fID, py, "py")
    print_prop(fID, pz, "pz")
    print_prop(fID, ρ, "ρ")
    print_prop(fID, Etot, "Etot")
    print_prop(fID, Ekin, "Ekin")
    print_prop(fID, Epot, "Epot")
    close(fID)
end

# Help function
# Function to write a single line
function print_prop(fID, x, sym)
    if typeof(x) == SingleDat
        spacestr = " "^(6-length(sym))
        if !isnan(x.val)
            @printf(fID,"%s:%s%.8e",sym,spacestr,x.val)
            if !isnan(x.std)
                @printf(fID," (%.8e",x.std)
                if !isnan(x.err)
                    @printf(fID,", %.8e",x.err)
                end
                @printf(fID,"),\n")
            else
                @printf(fID,",\n")
            end
        else
            @printf(fID,"%s:%s---,\n",sym,spacestr)
        end
    elseif typeof(x) == Array{SingleDat,1}
        for i = 1:length(x)
            print_prop(fID, x[i], string(sym,i))
        end
    elseif typeof(x) == Float64
        spacestr = " "^(6-length(sym))
        @printf(fID,"%s:%s%.8e",sym,spacestr,x)
        @printf(fID,",\n")
    end
end

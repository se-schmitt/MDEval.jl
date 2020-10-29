# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - EvalData
# Function to evaluate data from single MD run
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

# Function to ouput results of a single simulation / state point
function OutputResult(result::results_struct, folder::String)
    # Create Filepath
    path = string(folder,"/result.dat")

    # Write to file
    fID = open(path,"w")
    line1 = string("# Created by MD - Bulk Evaluation, Folder: ", folder)
    line2 = "# Format: val (std dev., std err.); [T]=K, [p]=MPa, [ρ]=g/ml, [E]=eV, [η]=Pa*s, [D]=m^2/s, [λ]=W/(m*K)"
    if (reduced_units) line2 = "# Format: val (std dev., std err.); reduced units" end
    header = string(line1,"\n",line2,"\n")
    print(fID,header)

    # Write data
    print_prop(fID, result.T, "T")
    print_prop(fID, result.p, "p")
    print_prop(fID, result.ρ, "ρ")
    print_prop(fID, result.Etot, "Etot")
    print_prop(fID, result.Ekin, "Ekin")
    print_prop(fID, result.Epot, "Epot")
    print_prop(fID, result.c, "c")
    print_prop(fID, result.η, "η")
    print_prop(fID, result.η_V, "η_V")
    print_prop(fID, result.D, "D")
    print_prop(fID, result.λ, "λ")
    close(fID)
end

# Function to ouput data from vle simulations (single simulation)
function OutputResult_VLE(dat::thermo_vle_dat, folder::String)
    # Create Filepath
    path = string(folder,"/result.dat")

    # Write to file
    fID = open(path,"w")
    line1 = string("# Created by MD - Bulk Evaluation, Folder: ", folder)
    line2 = "# Format: val (std dev., std err.); [T]=K, [p]=MPa, [ρ]=g/ml, [E]=eV"
    if (reduced_units) line2 = "# Format: val (std dev., std err.); reduced units" end
    header = string(line1,"\n",line2,"\n")
    print(fID,header)

    # Calculation of means, stds, errs
    T = single_dat(mean(dat.T),block_average(dat.T)[1],block_average(dat.T)[2])
    px = single_dat(mean(dat.px),block_average(dat.px)[1],block_average(dat.T)[2])
    py = single_dat(mean(dat.py),block_average(dat.py)[1],block_average(dat.T)[2])
    pz = single_dat(mean(dat.pz),block_average(dat.pz)[1],block_average(dat.T)[2])
    ρ = single_dat(mean(dat.ρ),block_average(dat.ρ)[1],block_average(dat.T)[2])
    Etot = single_dat(mean(dat.Etot),block_average(dat.Etot)[1],block_average(dat.T)[2])
    Ekin = single_dat(mean(dat.Ekin),block_average(dat.Ekin)[1],block_average(dat.T)[2])
    Epot = single_dat(mean(dat.Epot),block_average(dat.Epot)[1],block_average(dat.T)[2])

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
    if typeof(x) == single_dat
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
    end
end

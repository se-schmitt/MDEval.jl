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
    header = string(line1,"\n",line2,"\n")
    print(fID,header)

    # Write data
    print_prop(fID, result.T, "T")
    print_prop(fID, result.p, "p")
    print_prop(fID, result.ρ, "ρ")
    print_prop(fID, result.Etot, "Etot")
    print_prop(fID, result.Ekin, "Ekin")
    print_prop(fID, result.Epot, "Epot")
    print_prop(fID, result.η, "η")
    print_prop(fID, result.η_V, "η_V")
    print_prop(fID, result.D, "D")
    print_prop(fID, result.λ, "λ")
    close(fID)
end

# Help function
# Function to write a single line
function print_prop(fID, x, sym)
    if typeof(x) == single_dat
        spacestr = " "^(6-length(sym))
        if !isempty(x.val)
            if abs(ceil(log10(abs(x.val)))) > 3.5
                @printf(fID,"%s:%s%5.5e",sym,spacestr,x.val)
                if !isempty(x.std)
                    @printf(fID," (%5.5e",x.std)
                    if !isempty(x.err)
                        @printf(fID,", %5.5e",x.err)
                    end
                    @printf(fID,"),\n")
                else
                    @printf(fID,",\n")
                end
            else
                @printf(fID,"%s:%s%5.5f",sym,spacestr,x.val)
                if !isempty(x.std)
                    @printf(fID," (%5.5f",x.std)
                    if !isempty(x.err)
                        @printf(fID,", %5.5f",x.err)
                    end
                    @printf(fID,"),\n")
                else
                    @printf(fID,",\n")
                end
            end
        else
            @printf(fID,"%s:%s---,\n",sym,spacestr)
        end
    end
end

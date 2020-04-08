# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - EvalData
# Function to evaluate data from single MD run
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

# Function to ouput results of a single simulation / state point
function OutputResult(result, folder)
    # Create Filepath
    path = string(folder,"/result.dat")

    # Write to file
    fID = open(path,"w")
    line1 = string("# Created by MD - Bulk Evaluation, Folder: ", folder)
    line2 = "# Format: 'Quantity': val (std dev., std err.); [T]=K, [p]=MPa, [ρ]=g/ml, [E]=eV, [η]=Pa*s, [D]=m^2/s, [λ]=W/(m*K)"
    header = string(line1,"\n",line2,"\n")
    print(fID,header)

    if isempty(result.T.err) do_err = 0
    else                     do_err = 1 end

    # Write data
    if do_err == 0
        @printf(fID,"T:    %5.5f (%5.5f),\n", result.T.val, result.T.std)
        @printf(fID,"p:    %5.5f (%5.5f),\n", result.p.val, result.p.std)
        @printf(fID,"ρ:    %5.5f (%5.5f),\n", result.ρ.val, result.ρ.std)
        @printf(fID,"Etot: %5.5f (%5.5f),\n", result.Etot.val, result.Etot.std)
        @printf(fID,"Ekin: %5.5f (%5.5f),\n", result.Ekin.val, result.Ekin.std)
        @printf(fID,"Epot: %5.5f (%5.5f),\n", result.Epot.val, result.Epot.std)
        if !isempty(result.η.val)   @printf(fID,"η:    %5.5e,\n", result.η.val)
        else                        @printf(fID,"η:    ---,\n")     end
        if !isempty(result.η_V.val) @printf(fID,"η_V:  %5.5e,\n", result.η_V.val)
        else                        @printf(fID,"η_V:  ---,\n")     end
        if !isempty(result.D.val)   @printf(fID,"D:    %5.5e,\n", result.D.val)
        else                        @printf(fID,"D:    ---,\n")     end
        if !isempty(result.λ.val)   @printf(fID,"λ:    %5.5e,\n", result.λ.val)
        else                        @printf(fID,"λ:    ---,\n")     end
    elseif do_err == 1
        @printf(fID,"T:    %5.5f (%5.5f, %5.5f),\n", result.T.val, result.T.std, result.T.err)
        @printf(fID,"p:    %5.5f (%5.5f, %5.5f),\n", result.p.val, result.p.std, result.p.err)
        @printf(fID,"ρ:    %5.5f (%5.5f, %5.5f),\n", result.ρ.val, result.ρ.std, result.ρ.err)
        @printf(fID,"Etot: %5.5f (%5.5f, %5.5f),\n", result.Etot.val, result.Etot.std, result.Etot.err)
        @printf(fID,"Ekin: %5.5f (%5.5f, %5.5f),\n", result.Ekin.val, result.Ekin.std, result.Ekin.err)
        @printf(fID,"Epot: %5.5f (%5.5f, %5.5f),\n", result.Epot.val, result.Epot.std, result.Epot.err)
        if !isempty(result.η.val)
             @printf(fID,"η:    %5.5e (%5.5f, %5.5f),\n", result.η.val, result.η.std, result.η.err)
        else @printf(fID,"η:    ---,\n")     end
        if !isempty(result.η_V.val)
             @printf(fID,"η_V:  %5.5e (%5.5f, %5.5f),\n", result.η_V.val, result.η_V.std, result.η_V.err)
        else @printf(fID,"η_V:  ---,\n")     end
        if !isempty(result.D.val)
             @printf(fID,"D:    %5.5e (%5.5f, %5.5f),\n", result.D.val, result.D.std, result.D.err)
        else @printf(fID,"D:    ---,\n")     end
        if !isempty(result.λ.val)
             @printf(fID,"λ:    %5.5e (%5.5f, %5.5f),\n", result.λ.val, result.λ.std, result.λ.err)
        else @printf(fID,"λ:    ---,\n")     end
    end

    close(fID)
end

# Function to output the results of simulation at different state points
function OutputResult_multi(result, folder)
    # Create Filepath
    path = string(folder,"/result.dat")

    error("Not edited yet!")

    # Write to file
    fID = open(path,"w")
    units = "[T]=K, [p]=MPa, [ρ]=g/ml, [E]=eV, [η]=Pa*s, [D]=m^2/s, [λ]=W/(m*K)"
    header = string("Folder: ", folder,"\n",units,"\n---\n")
    print(fID,header)

    if isempty(result.T.err) do_err = 0
    else                     do_err = 1 end

    # Write data
    if do_err == 0
        @printf(fID,"T:    %5.5f (%5.5f),\n", result.T.val, result.T.std)
        @printf(fID,"p:    %5.5f (%5.5f),\n", result.p.val, result.p.std)
        @printf(fID,"ρ:    %5.5f (%5.5f),\n", result.ρ.val, result.ρ.std)
        @printf(fID,"Etot: %5.5f (%5.5f),\n", result.Etot.val, result.Etot.std)
        @printf(fID,"Ekin: %5.5f (%5.5f),\n", result.Ekin.val, result.Ekin.std)
        @printf(fID,"Epot: %5.5f (%5.5f),\n", result.Epot.val, result.Epot.std)
        if !isempty(result.η.val)   @printf(fID,"η:    %5.5e,\n", result.η.val)
        else                        @printf(fID,"η:    ---,\n")     end
        if !isempty(result.η_V.val) @printf(fID,"η_V:  %5.5e,\n", result.η_V.val)
        else                        @printf(fID,"η_V:  ---,\n")     end
        if !isempty(result.D.val)   @printf(fID,"D:    %5.5e,\n", result.D.val)
        else                        @printf(fID,"D:    ---,\n")     end
        if !isempty(result.λ.val)   @printf(fID,"λ:    %5.5e,\n", result.λ.val)
        else                        @printf(fID,"λ:    ---,\n")     end
    elseif do_err == 1
        @printf(fID,"T:    %5.5f (%5.5f, %5.5f),\n", result.T.val, result.T.std, result.T.err)
        @printf(fID,"p:    %5.5f (%5.5f, %5.5f),\n", result.p.val, result.p.std, result.p.err)
        @printf(fID,"ρ:    %5.5f (%5.5f, %5.5f),\n", result.ρ.val, result.ρ.std, result.ρ.err)
        @printf(fID,"Etot: %5.5f (%5.5f, %5.5f),\n", result.Etot.val, result.Etot.std, result.Etot.err)
        @printf(fID,"Ekin: %5.5f (%5.5f, %5.5f),\n", result.Ekin.val, result.Ekin.std, result.Ekin.err)
        @printf(fID,"Epot: %5.5f (%5.5f, %5.5f),\n", result.Epot.val, result.Epot.std, result.Epot.err)
        if !isempty(result.η.val)
             @printf(fID,"η:    %5.5e (%5.5f, %5.5f),\n", result.η.val, result.η.std, result.η.err)
        else @printf(fID,"η:    ---,\n")     end
        if !isempty(result.η_V.val)
             @printf(fID,"η_V:  %5.5e (%5.5f, %5.5f),\n", result.η_V.val, result.η_V.std, result.η_V.err)
        else @printf(fID,"η_V:  ---,\n")     end
        if !isempty(result.D.val)
             @printf(fID,"D:    %5.5e (%5.5f, %5.5f),\n", result.D.val, result.D.std, result.D.err)
        else @printf(fID,"D:    ---,\n")     end
        if !isempty(result.λ.val)
             @printf(fID,"λ:    %5.5e (%5.5f, %5.5f),\n", result.λ.val, result.λ.std, result.λ.err)
        else @printf(fID,"λ:    ---,\n")     end
    end

    close(fID)
end

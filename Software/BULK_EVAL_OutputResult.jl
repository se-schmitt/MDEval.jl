# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - EvalData
# Function to evaluate data from single MD run
# ---
# created by Sebastian Schmitt, 29.03.2020
# ------------------------------------------------------------------------------

# Main Function
function OutputResult(result, folder)
    # Create Filepath
    path = string(folder,"/result.dat")

    # Write to file
    fID = open(path,"w")
    units = "[T]=K, [p]=MPa, [ρ]=g/ml, [E]=eV, [η]=Pa*s, [D]=m^2/s, [λ]=W/(m*K)"
    header = string("Folder: ", folder,"\n",units,"\n---\n")
    print(fID,header)

    # Write data
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

    close(fID)
end

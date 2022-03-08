## NEMDShear.jl
# ------------------------------------------------------------------------------
# Evaluation Software for NEMD shear simulations (with sllod algorithm)
# Containing functions to evaluate output files of NEMD shear simulations to
# calculate viscosity
# ---
# created by Sebastian Schmitt, 16.07.2021

function EvalNEMDShear(subfolder,inpar)
    moltype, dt, natoms, molmass = load_info(subfolder)

    # Initialization of info structure
    info = info_struct( subfolder,          # info.folder
                        inpar.ensemble,     # info.ensemble
                        inpar.n_equ,        # info.n_equ
                        moltype,            # info.moltype
                        dt,                 # info.dt
                        natoms,             # info.natoms
                        molmass,            # info.molmass
                        inpar.n_blocks,     # info.n_blocks
                        inpar.N_bin,        # info.N_bin
                        inpar.r_cut)        # info.r_cut

    # Average Thermodynamic Properties
    T, p, ρ, Etot, Ekin, Epot, c, pyz = ave_thermo(info; is_nemd=true)

    # Load data
    ts_add = 0
    data = []
    for f in list[startswith.(list,"vy_profile")]
        data, no_chunks, n_steps = read_profile1D("$(info.folder)/$f",prof,ts_add)
        ts_add = prof.timestep[end]
    end
    what = data.timestep .>= info.n_equ
    vy_mean = return_vy_mean(data,what,n_steps,no_chunks)
    x = vy_mean[:,2]
    y = vy_mean[:,4]
    a,b = hcat(fill!(similar(x), 1), x) \ y
    x_0 = ones(size(x,1))

    y_predict=a*x_0+b*x
    # Sum of Squares Total, Sum of Squares Regression
    SST = sum((y .- mean(y)).^2)
    SSR = sum((y_predict .- mean(y)).^2)
    r_square = SSR / SST
    r_squared = single_dat(r_square,0,0)


    timestep_equend = data.timestep[what]
    filename2 = "$(info.folder)/in.Bulk_NVT_sllod"
    N_mol = read_inputfile(filename2)

    L = ((N_mol*info.molmass*1e24)/(ρ.val*6.02214076e23))^(1/3)
    s_rate = 10^12*b/L #10^15 for ReaxFF 10^12 for others
    dat = load_thermo(info, is_nemd="shear")
    if (reduced_units)      factor_p = 1
    elseif !(reduced_units) factor_p = 0.1 end

    η_vec = -(dat.pyz.*factor_p)*1e6/s_rate

    η = single_dat(mean(η_vec), block_average(η_vec; M_block=100)[1], block_average(η_vec; M_block=100)[2])
    figure()
    xlabel(L"x")
    ylabel(L"vy")
    plot(x,y, color="gray", linestyle=":", linewidth=0.1)
    title("Fit",fontsize=10)
    legend(loc="upper left")
    tight_layout()
    filename3 = "$(info.folder)/vy_plot.pdf"
    savefig(filename3)

    # Output results
    res = results_struct_nemd(T, p, ρ, x, Etot, Ekin, Epot, pyz, η, s_rate, r_squared, [])
    OutputResultNEMD(res, info.folder)
end




## Subfunctions
function return_vy_mean(data,what,n_steps,no_chunks)

    id_chunk_mat = ones(1,n_steps)
    x_mat = ones(1,n_steps)
    Ncount_mat = ones(1,n_steps)
    vy_mat = ones(1,n_steps)
    vy_mean = ones(no_chunks,6)


    #loop to generate vy_mean matrix for n_equ:end
    for i = 1 : no_chunks
        for j = 1 : n_steps
            id_chunk_mat[j] = data.id_chunk[j,i]
            x_mat[j] = data.x[j,i]
            Ncount_mat[j] = data.Ncount[j,i]
            vy_mat[j] = data.vy[j,i]
        end
        id_chunk_vec = single_dat(mean(id_chunk_mat[what]), block_average(id_chunk_mat[what])[1], block_average(id_chunk_mat[what])[2])
        x_vec = single_dat(mean(x_mat[what]), block_average(x_mat[what])[1], block_average(x_mat[what])[2])
        Ncount_vec = single_dat(mean(Ncount_mat[what]), block_average(Ncount_mat[what])[1], block_average(Ncount_mat[what])[2])
        vy_vec = single_dat(mean(vy_mat[what]), block_average(vy_mat[what])[1], block_average(vy_mat[what])[2])

        #stores mean values for n_equ:end for col1: id_chunk, col2: x, col3:N_count, col4-6 vy
        vy_mean[i,1] = round(id_chunk_vec.val)
        vy_mean[i,2] = round(x_vec.val,digits=4)
        vy_mean[i,3] = round(Ncount_vec.val,digits=0)
        vy_mean[i,4] = round(vy_vec.val,digits=4)
        vy_mean[i,5] = round(vy_vec.std,digits=4)
        vy_mean[i,6] = round(vy_vec.err,digits=4)
    end
    return vy_mean
end

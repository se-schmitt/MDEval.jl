## StokesEinstein.jl
# ------------------------------------------------------------------------------
# Evaluation Software for MD Bulk Simulations - StokesEinstein
# Containing functions to apply Stokes-Einstein method to self difussion
# coefficient
# ---
# created by Sebastian Schmitt, 15.03.2021
# ------------------------------------------------------------------------------

## Self diffusion coefficient --------------------------------------------------
# Evaluate Atoms Positions to calculate Self Diffusion Coefficient
function calc_selfdiffusion(info::info_struct, dat::Array{Any,1})
    if (!isempty(dat))
        N_moltype = maximum(dat[1].moltype)
        D = Array{single_dat,1}(undef,N_moltype)

        # Convert atom to molecule coordinates
        mol = atom2mol(dat)

        # Initialize figure and msd_t_all variable
        colors = ["b","g","r","c","m","y"]
        figure()
        msd_t_all = Array{Float64,2}(undef,0,N_moltype)

        for i = 1:N_moltype
            # Get all molecules of type 'moltype'
            mol_i = get_mol_by_type(mol,i)

            # Calculation of Mean Square Displacement
            msd_t, t = calc_msd(mol_i,info)

            # Segmenting the msd and calculation of log-log slope
            kmax_eval = 0.6
            δ_tol = 0.075
            δ = 0
            N_eval = ceil(Int64,kmax_eval*length(t))
            L_eval = 50
            N = N_eval
            while δ < δ_tol && N > L_eval
                N = N-L_eval
                what = Array(N+1:N+L_eval)
                beta = hcat(ones(L_eval),log10.(t[what]))\log10.(msd_t[what])
                δ = abs(beta[2]-1)
            end

            what_eval = Int64.(N+1:N_eval)
            if length(what_eval) < 150
                what_eval = Int64.(ceil(0.25*length(t)):ceil(0.75*length(t)))
            end

            # Calculation of D by determination of slope of msd
            beta = hcat(ones(length(what_eval)),t[what_eval]) \ msd_t[what_eval]
            factor_unit = 1e-8
            if (reduced_units) factor_unit = 1 end

            D[i] = single_dat(beta[2]/6*factor_unit,NaN,NaN)

            # Plot MSD of molecule i
            t1 = round(Int64,t[what_eval[1]])
            t2 = round(Int64,t[what_eval[end]])

            plot(t, msd_t, linewidth=1, color=colors[i], label=string("Mol. ",i," (",t1," ps - ",t2," ps (",length(what_eval)," steps))"))
            plot(t,beta[1] .+ beta[2].*t, linestyle="-", linewidth=1, color=colors[i])
            plot([t1,t1],[0.1,msd_t[what_eval[1]]], color="black", linestyle=":", linewidth=1)
            plot([t2,t2],[0.1,msd_t[what_eval[end]]], color="black", linestyle=":", linewidth=1)

            # Save MSD in array to save dlm file
            if i == 1
                msd_t_all = t
            end
            msd_t_all = hcat(msd_t_all, msd_t)
        end

        # Save msd as data file
        line1 = string("# Created by MD - Bulk Evaluation, Folder: ", info.folder)
        line2 = string("# t[ps]"," msd[Å]"^N_moltype)
        header = string(line1,"\n",line2)
        file = string(info.folder,"/msd.dat")
        fID = open(file,"w"); println(fID,header)
        writedlm(fID, msd_t_all, " ")
        close(fID)

        # Formatting figure
        if !(reduced_units) xlabel(L"t / ps") elseif reduced_units xlabel(L"t*") end
        if !(reduced_units) ylabel(L"MSD / Å$^2$") elseif reduced_units ylabel(L"MSD*") end
        title("Mean square displacement (MSD)")
        legend()
        tight_layout()
        savefig(string(info.folder,"/fig_msd(t).pdf"))
        close()
    else
        D = single_dat(NaN,NaN,NaN)
    end
    return D
end

# Function to calculate mean square displacement
function calc_msd(dat::Array{Any,1},info)
    n = length(dat[1].molid)    # No. molecules
    N = length(dat)             # No. timesteps

    t = Float64[]
    for tmp in dat append!(t,tmp.t) end
    x = fill(Float64[],n)
    y = fill(Float64[],n)
    z = fill(Float64[],n)
    i = 0
    for imol in dat[1].molid
        i += 1
        xmol = Float64[]
        ymol = Float64[]
        zmol = Float64[]
        for tmp in dat
            what = tmp.molid .== imol
            append!(xmol,tmp.x[what])
            append!(ymol,tmp.y[what])
            append!(zmol,tmp.z[what])
        end
        x[i] = xmol; y[i] = ymol; z[i] = zmol
    end

    msd_t = zeros(N)

    for i = 1:n
        r2 = x[i] .^ 2 + y[i] .^ 2 + z[i] .^ 2

        s_AB = sum(pmap(x -> calc_acf(x,"fft"),[x[i],y[i],z[i]]))

        msd_i = Array{Float64,1}(undef,N)
        sumsq = []
        for im = 1:N
            m = im-1
            if (im == 1) sumsq = 2*sum(r2)
            else         sumsq = sumsq - r2[im-1] - r2[N-m] end
            msd_i[im] = sumsq / (N - m) - 2*s_AB[im]
        end

        msd_t += msd_i ./n
    end
    msd_t = abs.(msd_t)

    return msd_t, t
end

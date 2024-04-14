# Struct to save evaluation settings
mutable struct Opts
    mode::String
    folder::String
    ensemble::String
    n_equ::Union{Int64,Missing}
    do_single::Union{Bool,Missing}
    do_state::Union{Bool,Missing}
    n_boot::Union{Int64,Missing}
    cutcrit::Union{Float64,Missing}
    do_transport::Union{Bool,Missing}
    corr_length::Union{Int64,Missing}
    span_corr_fun::Union{Int64,Missing}
    n_blocks::Union{Int64,Missing}
    n_every::Union{Int64,Missing}
    debug_mode::Union{Int64,Missing}
    acf_calc_mode::Union{String,Missing}
    do_structure::Union{Bool,Missing}
    N_bin::Union{Int64,Missing}
    r_cut::Union{Float64,Missing}
    k_L_thermo::Union{Float64,Missing}
    Opts() = new("","","",repeat([missing],16)...)
end

# Function to convert input to Opts structure
function get_opts(mode, folder, kw)
    # Initialization of input variables
    opts = Opts()

    # Set units
    if :units in keys(kw) && kw.units == "reduced"
        reduced_units = true
    end

    # Set options 'general'
    opts.mode = String(mode)
    opts.folder = folder
    opts.ensemble = :ensemble in keys(kw) ? kw.ensemble : ""
    opts.n_equ = :timesteps_equ in keys(kw) ? kw.timesteps_equ : 0
    opts.n_blocks = :n_blocks in keys(kw) ? kw.n_blocks : 0
    opts.do_single = :do_single in keys(kw) ? kw.do_single : true
    opts.do_state = :do_state in keys(kw) ? kw.do_state : true
    opts.debug_mode = :debug_mode in keys(kw) ? kw.debug_mode : false
    opts.do_structure = :do_structure in keys(kw) ? kw.do_structure : false
    if opts.do_structure
        opts.N_bin = :n_bin in keys(kw) ? kw.n_bin : 100
        opts.r_cut = :r_cut in keys(kw) ? kw.r_cut : 10.0
    end

    # Set options 'single_run'
    if mode == :single_run
        opts.do_transport = :do_transport in keys(kw) ? kw.do_transport : true
        if opts.do_transport
            opts.corr_length = :corr_length in keys(kw) ? kw.corr_length : throw_missing(:corr_length)
            opts.span_corr_fun = :span_corr_fun in keys(kw) ? kw.span_corr_fun : throw_missing(:span_corr_fun)
        end
        opts.n_every = :n_every in keys(kw) ? kw.n_every : 1
        opts.acf_calc_mode = :acf_calc_mode in keys(kw) ? kw.acf_calc_mode : "autocov"
    end

    # Set options 'tdm'
    if mode == :tdm
        opts.n_boot = :n_boot in keys(kw) ? kw.n_boot : throw_missing(:n_boot)
        opts.cutcrit = :cutcrit in keys(kw) ? kw.cutcrit : 0.4
        opts.acf_calc_mode = :acf_calc_mode in keys(kw) ? kw.acf_calc_mode : "fft"
    end

    # Set options 'nemd'
    if mode == :nemd_shear
        opts.ensemble = "NVT"
    end
    if mode == :nemd_heat
        opts.ensemble = "NVE"
        opts.k_L_thermo = :k_L_thermo in keys(kw) ? kw.k_L_thermo : 0.1
        opts.n_blocks = :n_blocks in keys(kw) ? kw.n_blocks : 20
    end

    return opts
end

function throw_missing(name)
    error("Parameter ':$name' is missing!")
end
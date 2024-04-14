# Example: Singe Run

using MDEval

mode = :single_run
folder = "./EMD_LJ"
keywords = (;
    ensemble        =   "NVT",
    do_transport    =   true,
    corr_length     =   20000,
    span_corr_fun   =   500,
    units           =   "reduced",   
)

mdeval(mode,folder,keywords)

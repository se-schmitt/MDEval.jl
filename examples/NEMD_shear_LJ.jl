# Example: NEMD shear

using MDEval

mode = :nemd_shear
folder = "./NEMD_shear_LJ/"
keywords = (;
    do_state        =   false,
    units           =   "reduced",   
)

mdeval(mode,folder,keywords)

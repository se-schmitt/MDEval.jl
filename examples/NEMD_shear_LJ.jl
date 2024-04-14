# Example: NEMD shear

using MDEval

mode = :nemd_shear
folder = "./NEMD_shear_LJ/"
keywords = (;
    units           =   "reduced",   
)

mdeval(mode,folder,keywords)

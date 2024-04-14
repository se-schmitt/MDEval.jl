# Example: NEMD heat

using MDEval

mode = :nemd_heat
folder = "./NEMD_heat_LJ/"
keywords = (;
    units           =   "reduced",   
)

mdeval(mode,folder,keywords)
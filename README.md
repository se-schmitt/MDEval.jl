# Documentation - MD_Evaluation

Program to evaluate MD output files of bulk simulations
Notes:

- usable for LAMMPS simulations in metal units
- required files provided by specific LAMMPS scripts

## 1. Overview

Program evaluates four different simulation types:

- Evaluation of **single run simulations** for transport properties (**single_run**)
- **Time decomposition method (TDM)** for transport properties (**tdm**)
- Evaluation of **VLE simulations** (two phase simualtions) (**vle**)
- Evaluation of **NEMD shear simulations** (**nemd-shear**)

### Single run simulations

**Required Folder structure:**

- *SIM_1* [contains all output files of a EMD simulation of one state point]
- ...

### TDM

**Required Folder structure:**

- *STATE_1* [contains all simulations of one state point]
  - *Sim_001*
  - *Sim_002*
  - ...
- ...

### VLE simulations

**Required Folder structure:**

- *FOLDER_1* [contains all simulations done for one substance]
  - *Sim_T1* [VLE simulation at temperature T1]
  - *Sim_T2* [VLE simulation at temperature T2]
  - ...
- ...

### NEMD shear simulations

**Required Folder structure:**

- *STATE_1* [contains all simulations of one state point]
  - *Sim_S1* [VLE simulation with shear rate S1]
  - *Sim_S2* [VLE simulation with shear rate S2]
  - ...
- ...

## 2. Input File

**Selecting a input file:**

Two options:

- There has to be a file *INPUT.txt*  in the MD_Evaluation folder **or**
- the path of the file can be passed as argument when running the master function *MD_EVAL_Master.jl* (julia *MD_EVAL_Master.jl* [*name of the input file*])

**General structure of the input file:**

`name_of_keyword = value_for_parameter`

### Keywords:

| Name              | Type {*standard value*}                                         | Modes                | Description                                                                                                                                                                       |
| ----------------- | --------------------------------------------------------------- | -------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **mode**          | string [*single_run*, *tdm*, *vle*, *nemd-shear*]               | all                  | defines the mode of the simulations/evaluation                                                                                                                                    |
| **folder**        | string                                                          | all                  | path to main folder containing all simulation data (see chapter 1)                                                                                                                |
| **ensemble**      | string [*NVT*, *NVE*, *NpT*]                                    | single_run, tdm      | ensemble to evaluate                                                                                                                                                              |
| **timesteps_equ** | interger (≥ 0)                                                  | single_run, tdm, vle | number of timesteps to ignore at the start of each simulation                                                                                                                     |
| **do_single**     | integer [0,1]                                                   | tdm                  | 1 - evaluate single folders, 0 - single folders already evaluated                                                                                                                 |
| **do_state**      | integer [0,1]                                                   | tdm                  | 1 - evaluate complete thermodynamic state (main folder)                                                                                                                           |
| **n_boot**        | integer (≥ 0)                                                   | tdm                  | number of bootstrapping repetitions                                                                                                                                               |
| **cutcrit**       | float (≥ 0) {*0.4*}                                             | tdm                  | cut criteria for tdm method                                                                                                                                                       |
| **do_transport**  | integer [0,1]                                                   | single_run           | 1 - do evaluation of transport properties, 0 - skip evaluation of transport properties                                                                                            |
| **corr_length**   | integer (≥ 0)                                                   | single_run           | length (timesteps) of correlation function                                                                                                                                        |
| **span_corr_fun** | integer (≥ 0)                                                   | single_run           | timesteps between single correlation functions                                                                                                                                    |
| **n_blocks**      | integer (≥ 0)                                                   | single_run           | number of blocks for static properties                                                                                                                                            |
| **n_every**       | integer (≥ 1)                                                   | single_run           | skip n_every timesteps when calculating acf (useful for slowly converging states, e.g. ideal gas)                                                                                 |
| **debug_mode**    | integer [0,1] {*0*}                                             | single_run           | if debug mode is enables, errors are thrown directly                                                                                                                              |
| **acf_calc_mode** | string [*autocov*, *fft*] {single_run → *autocov*, tdm → *fft*} | single_run, tdm      | mode for acf calculation (*autocov*: full acf by Julia *autocov* command (can be slowly for long signals), *fft*: acf calculation by FFT (fast, but inaccurate for long signals)) |
| **do_structure**  | integer [0,1] {*0*}                                             | all                  | 1 - do structure evaluation, 0 - skip structure evaluation                                                                                                                        |
| **n_bin**         | integer (≥ 0) {*100*}                                           | all                  | number of bins for rdf calculation                                                                                                                                                |
| **r_cut**         | float (unit: Å) {*10 Å*}                                        | all                  | cut-off radius for rdf calculation                                                                                                                                                |
| **units**         | string [*real*, *reduced*] {*real*}                             | all                  | units of simulation (real: LAMMPS SI units, reduced: reduced by LJ parameters)                                                                                                    |
| **k_L_thermo**    | float ($0\leq k \leq 0.5$) {0.25}                               | nemd-heat            | reduced length of thermostats ($L_x^*=1$)                                                                                                                                         |

Example *single*:

```
# Mode
mode          =   single_run

# Folder
folder        =   C:/path2simulations/sim_1
folder        =   C:/path2simulations/sim_*

# Ensemble and equilibration
ensemble      =   NVT
timesteps_EQU =   0

# Settings for acf
corr_length   =   100000
span_corr_fun =   20000
```

Example *TDM*:

```
# Mode
mode          =   tdm

# Folder
folder        =   C:/path2simulations/sim_1
folder        =   C:/path2simulations/sim_*

# Ensemble and equilibration
ensemble      =   NVT
timesteps_EQU =   0

# TDM Settings
DO_single     =   1
DO_state      =   1
N_boot        =   100
```

## Output

### Units

| Property                   | Symbol | Unit    |
| -------------------------- | ------ | ------- |
| mass                       | *m*    | g/mol   |
| time                       | *t*    | ps      |
| temperature                | *T*    | K       |
| pressure                   | *p*    | MPa     |
| energy                     | *E*    | eV      |
| viscosity                  | *η*    | Pa*s    |
| self-diffusion coefficient | *D*    | m²/s    |
| thermal conductivity       | *λ*    | W/(m*K) |

## Issues

- [ ] NEMD parts  
  - [ ] provide $x$ (mole fraction) 
  - [ ] improve error calculation 
- [ ] info File: reading the file should not depend on line numbers

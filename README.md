# Usage


# To-Dos

- [ ] add command line argument to choose input file (optional)


# Documentation - MD_Evaluation

Program to evaluate MD output files of bulk simulations
Notes:
- usable for LAMMPS simulations in metal units
- required files provided by specific LAMMPS scripts

## 1. Overview

Programm evaluates three different simulation types:
- Evaluation of **single run simulations** for transport properties (**single_run**)
- **Time decomposition method (TDM)** for transport properties (**tdm**)
- Evaluation of **VLE simulations** (two phase simualtions) (**vle**)

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

## 2. Input File

**Name**: 'INPUT.txt' (has to be in same folder as 'Master' file)

**General structure**:
  - first line: '#' + 'keyword'
  - second line: value assigned to keyword

### Keywords:

| Name              | Type {*standard value*}           | Modes           | Description |
| ----------------- | --------------------------------- | --------------- | ----------- |
| **mode**          | string [*single_run*, *tdm*, *vle*] | all             | defines the mode of the simulations/evaluation |
| **folder**        | string                            | all             | path to main folder containing all simulation data (see chapter 1) |
| **ensemble**      | string [*NVT*, *NVE*, *NpT*]        | single_run, tdm | ensemble to evaluate |
| **timesteps_EQU** | interger (≥ 0)                    | all         | number of timesteps to ignore at the start of each simulation |
| **DO_single**     | integer [0,1]                     | tdm         | 1 - evaluate single folders, 0 - single folders already evaluated |
| **DO_state**      | integer [0,1]                     | tdm         | 1 - evaluate complete thermodynamic state (main folder) |
| **N_boot**        | integer (≥ 0)                     | tdm         | number of bootstrapping repetitions |
| **DO_transport** | integer [0,1] | single_run | 1 - do evaluation of transport properties, 0 - skip evaluation of transport properties |
| **corr_length**   | integer (≥ 0)                     | single_run  | length (timesteps) of correlation function |
| **span_corr_fun** | integer (≥ 0)                     | single_run  | timesteps between single correlation functions |
| **n_every**       | integer (≥ 1)                     | single_run  | skip n_every timesteps when calculating acf (useful for slowly converging states, e.g. ideal gas) |
| **DO_structure**  | integer [0,1] {*0*}               | all         | 1 - do structure evaluation, 0 - skip structure evaluation |
| **N_bin**         | integer (≥ 0) {*100*}             | all         | number of bins for rdf calculation |
| **r_cut**         | float (unit: Å) {*10 Å*}          | all         | cut-off radius for rdf calculation |
| **units**         | string [*real*, *reduced*] {*real*} | all       | units of simulation (real: LAMMPS SI units, reduced: reduced by LJ parameters) |

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

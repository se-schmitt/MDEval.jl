# MDEval.jl

Package to evaluate molecular dynamics (MD) simulations.

> [!WARNING]  
> Work in progress. Currently not working as documented here.

## 1. Overview

### 1.1 Installation

To register the module locally, type 
```julia
Pkg> add https://github.com/se-schmitt/MDEval.jl
```
in package mode (type `]` in REPL to enter Pkg mode).

Then, the module can be loaded by
```julia
using MDEval
```

### 1.2 Running an Evaluation

```julia
using MDEval

mode = :single_run
folder = "your/simulation/folder"
keywords = (;
    ensemble        =   "NVT",
    timesteps_equ   =   1e5,
    do_transport    =   true,
    corr_length     =   100000
    span_corr_fun   =   20000
)

MDeval(mode,folder,keywords)
```

### 1.3. Evaluation Modes

Program evaluates four different simulation types:

- Evaluation of single run simulations for thermodynamic properties (including transport properties) (`:single_run`)
- [Time decomposition method (TDM)](https://www.doi.org/10.1021/acs.jctc.5b00351) for transport properties (`:tdm`)
- Evaluation of VLE simulations (direct two phase simualtions) (`:vle`)
- Evaluation of NEMD shear simulations to determine viscosity (`:nemd_shear`)
- Evaluation of NEMD heat transfer simulation to determine thermal conductivity (`:nemd_heat`)

## 2. Manual

### 2.1 Function `MDEval`

`MDEval(mode::Symbol,folder::String,keywords::NamedTuple)`

Input arguments:
- `mode`: defines the mode of the simulations/evaluation. Possible values: `:single_run`, `:tdm,:vle`, `:nemd_shear`, or `:nemd_heat`
- `folder`: path to main folder containing all simulation data (see [Section 2.2](#22-folder-structure))
- `keywords`: options to define the evaluation parameters (see [Section 2.1](#21-keywords))

### 2.1 Keywords

| Name              | Possible values [default value]                                         | Modes                | Description                                                                                                                                                                       |
| ----------------- | --------------------------------------------------------------- | -------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
| `ensemble`      | `NVT`, `NVE`, `NpT`                                     | `:single_run`, `:tdm`, `:nemd_shear`     | ensemble to evaluate                                                                                                                                                              |
| `timesteps_equ` | interger ≥ 0 [`0`]                                                  | `:single_run`, `:tdm`, `:vle` | number of timesteps to ignore at the start of each simulation                                                                                                                     |
| `do_single`     | `true`, `false` [`true`]                                                  | `:tdm`, `:nemd_shear`      | 1 - evaluate single folders, 0 - single folders already evaluated                                                                                                                 |
| `do_state`      | `true`, `false` [`true`]                                                  | `:tdm`, `:nemd_shear`      | 1 - evaluate complete thermodynamic state (main folder)                                                                                                                           |
| `n_boot`        | integer ≥ 0 [`0`]                                                   | `:tdm`                  | number of bootstrapping repetitions                                                                                                                                               |
| `cutcrit`       | float ≥ 0 [`0.4`]                                             | `:tdm`                  | cut criteria for `:tdm` method                                                                                                                                                       |
| `do_transport`  | `true`, `false`                                                    | `:single_run`           | 1 - do evaluation of transport properties, 0 - skip evaluation of transport properties                                                                                            |
| `corr_length`   | integer ≥ 0                                                   | `:single_run`           | length (timesteps) of correlation function                                                                                                                                        |
| `span_corr_fun` | integer ≥ 0                                                   | `:single_run`           | timesteps between single correlation functions                                                                                                                                    |
| `n_blocks`      | integer ≥ 0                                                   | `:single_run`           | number of blocks for static properties                                                                                                                                            |
| `n_every`       | integer ≥ 1                                                   | `:single_run`           | skip n_every timesteps when calculating acf (useful for slowly converging states, e.g. ideal gas)                                                                                 |
| `debug_mode`    | `true`, `false`                                             | `:single_run`           | if debug mode is enables, errors are thrown directly                                                                                                                              |
| `acf_calc_mode` | `autocov`, `fft` [`:single_run`→`autocov`, `:tdm`→`fft`] | `:single_run`, `:tdm`      | mode for acf calculation (`autocov`: full acf by Julia `autocov` command (can be slowly for long signals), `fft`: acf calculation by FFT (fast, but inaccurate for long times)) |
| `do_structure`  | `true`, `false` [`false`]                                             | all                  | 1 - do structure evaluation, 0 - skip structure evaluation                                                                                                                        |
| `n_bin`         | integer ≥ 0 [`100`]                                           | all                  | number of bins for rdf calculation                                                                                                                                                |
| `r_cut`         | float [`10.0`]                                        | all                  | cut-off radius for rdf (unit: Å)calculation                                                                                                                                                |
| `units`         | `metal`, `reduced` [`real`]                             | all                  | units of simulation (real: LAMMPS SI units, reduced: reduced by LJ parameters)                                                                                                    |
| `k_L_thermo`    | float [`0.1`]                               | `:nemd_heat`            | reduced length of thermostats ($L_x^*=1$)                                                                                                                                         |

### 2.2 Folder structure

#### 2.2.1 Single run simulations

```
├── SIM_1* [contains all output files of a EMD simulation of one state point]
:
```

#### 2.2.2 TDM

```
├── STATE_1 [contains all replica simulations of one state point]
│   ├── Sim_001
│   ├── Sim_002
:   :
```

#### 2.2.3 VLE simulations

```
├── FOLDER_1 [contains all simulations for one substance]
│   ├── Sim_T1 [VLE simulation at temperature T1]
│   ├── Sim_T2 [VLE simulation at temperature T2]
:   :
```

#### 2.2.4 NEMD shear simulations

```
├── STATE_1 [contains all simulations of one state point]
│   ├── Sim_S1 [VLE simulation with shear rate S1]
│   ├── Sim_S2 [VLE simulation with shear rate S2]
:   : 
```

## 3. Output

## 4. Units

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

## 5. Parallel execution

## Issues

- [ ] make it a module
- [ ] NEMD parts  
  - [ ] provide $x$ (mole fraction)
  - [ ] improve error calculation
- [ ] info File: reading the file should not depend on line numbers
- [ ] remove ensemble from filenames that do not need it (e.g. nemd-shear)

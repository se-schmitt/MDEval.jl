# MDEval.jl

Package to evaluate molecular dynamics (MD) simulations.

## 1. Installation

To register the module locally, type 
```julia
Pkg> add https://github.com/se-schmitt/MDEval.jl
```
in package mode (type `]` in REPL to enter Pkg mode).

Then, the module can be loaded by
```julia
using MDEval
```

## 2. Manual

### 2.1 Function `mdeval`

`mdeval(mode::Symbol, folder::String, keywords::NamedTuple)`

Input arguments:
- `mode`: defines the mode of the simulations/evaluation. Possible values: `:single_run`, `:tdm,:vle`, `:nemd_shear`, or `:nemd_heat` (see [Section 2.1](#21-evaluation-modes))
- `folder`: path to main folder containing all simulation data (see [Section 2.3](#23-folder-structure))
- `keywords`: options to define the evaluation parameters (see [Section 2.2](#22-keywords))

Examples are provided in the `examples` folder.

### 2.1 Evaluation Modes

Different simulation types can be evaluated:

- Evaluation of single run simulations for thermodynamic properties (including transport properties) (`:single_run`)
- [Time decomposition method (TDM)](https://www.doi.org/10.1021/acs.jctc.5b00351) for transport properties (`:tdm`)
- Evaluation of VLE simulations (direct two phase simualtions) (`:vle`)
- Evaluation of NEMD shear simulations to determine viscosity (`:nemd_shear`)
- Evaluation of NEMD heat transfer simulation to determine thermal conductivity (`:nemd_heat`)

### 2.2 Keywords

| Name            | Possible values [default value]    | Modes                                | Description                                                                                                                                                                     |
|-----------------|------------------------------------|--------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `ensemble`      | `NVT`, `NVE`, `NpT`                | `:single_run`, `:tdm`, `:nemd_shear` | ensemble to evaluate                                                                                                                                                            |
| `timesteps_equ` | interger ≥ 0 [`0`]                 | `:single_run`, `:tdm`, `:vle`        | number of timesteps to ignore at the start of each simulation                                                                                                                   |
| `do_single`     | `true`, `false` [`true`]           | `:tdm`, `:nemd_shear`                | evaluate single folders                                                                                                                                                         |
| `do_state`      | `true`, `false` [`true`]           | `:tdm`, `:nemd_shear`                | evaluate complete thermodynamic state (main folder)                                                                                                                             |
| `n_boot`        | integer ≥ 0 [`0`]                  | `:tdm`                               | number of bootstrapping repetitions                                                                                                                                             |
| `cutcrit`       | float ≥ 0 [`0.4`]                  | `:tdm`                               | cut criteria for `:tdm` method                                                                                                                                                  |
| `do_transport`  | `true`, `false`                    | `:single_run`                        | evaluation of transport properties                                                                                                                                              |
| `corr_length`   | integer ≥ 0                        | `:single_run`                        | length (timesteps) of correlation function                                                                                                                                      |
| `span_corr_fun` | integer ≥ 0                        | `:single_run`                        | timesteps between single correlation functions                                                                                                                                  |
| `n_blocks`      | integer ≥ 0                        | `:single_run`                        | number of blocks for static properties                                                                                                                                          |
| `n_every`       | integer ≥ 1                        | `:single_run`                        | skip n_every timesteps when calculating acf (useful for slowly converging states, e.g. ideal gas)                                                                               |
| `debug_mode`    | `true`, `false`                    | `:single_run`                        | if debug mode is enables, errors are thrown directly                                                                                                                            |
| `acf_calc_mode` | `autocov`, `fft` [`autocov`/`fft`] | `:single_run`, `:tdm`                | mode for acf calculation (`autocov`: full acf by Julia `autocov` command (can be slowly for long signals), `fft`: acf calculation by FFT (fast, but inaccurate for long times)) |
| `do_structure`  | `true`, `false` [`false`]          | all                                  | structure evaluation                                                                                                                                                            |
| `n_bin`         | integer ≥ 0 [`100`]                | all                                  | number of bins for rdf calculation                                                                                                                                              |
| `r_cut`         | float [`10.0`]                     | all                                  | cut-off radius for rdf (unit: Å)calculation                                                                                                                                     |
| `units`         | `metal`, `reduced` [`real`]        | all                                  | units of simulation (real: LAMMPS SI units, reduced: reduced by LJ parameters)                                                                                                  |
| `k_L_thermo`    | float [`0.1`]                      | `:nemd_heat`                         | reduced length of thermostats (as fraction of total box length)                                                                                                                 |

### 2.3 Folder structure

- Single run/ NEMD shear/ NEMD heat simulations
```
├── SIM_1* [contains all output files of a EMD simulation of one state point]
:
```

- TDM
```
├── STATE_1 [contains all replica simulations of one state point]
│   ├── Sim_001
│   ├── Sim_002
:   :
```

- VLE simulations
```
├── FOLDER_1 [contains all simulations for one substance]
│   ├── Sim_T1 [VLE simulation at temperature T1]
│   ├── Sim_T2 [VLE simulation at temperature T2]
:   :
```

- NEMD shear simulations
```
├── STATE_1 [contains all simulations of one state point]
│   ├── Sim_S1 [VLE simulation with shear rate S1]
│   ├── Sim_S2 [VLE simulation with shear rate S2]
:   : 
```

### 2.4 Units

Either `reduced` or `real` units can be used. The given unit defines the input as well as the output units. `reduced` units refer to the Lennard-Jones units system. The input units for `real` units are the LAMMPS [metal units](https://www.afs.enea.it/software/lammps/doc17/html/units.html). The units of the output are defined in the following table.

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

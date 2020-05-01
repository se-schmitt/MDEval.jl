
# MD - Bulk Evaluation

Program to evaluate MD output files of bulk simulations
Notes:
- usable for LAMMPS simulations in metal units
- required files provided by LAMMPS scipt "in.output"

## Manual

__Required Folder structure:__
- "STATE_1"
  - "Sim_001"
  - "Sim_002"
  - ...
- "STATE_1"
  - "Sim_001"
  - "Sim_002"
  - ...
- ...

## Input

__Input Variables:__

| Name      | Possible Values             | Description                                                           |
| --------- | ----------------------------| --------------------------------------------------------------------- |
| folders   | ["STATE_1", "STATE_2", ...] | Array of strings with paths to simulation folders                     |
| ensemble  | "NVT" / "NVE" / "NpT"       | String representing ensemble to evaluate                              |
| n_equ     | 10000                       | Integer with number of equilibration steps (are skipped at beginning) |
| do_eval   | 1 / 0                       | 1: Evaluation of single simulations; 0: Simulations already evaluated |
| do_state  | 1 / 0                       | 1: Calculation of state variables on basis of single Simulations      |
| nboot     | 200                         | Number of bootstrapping loop                                          |
| no_procs  | 8                           | Number of processors to use for evaluation                            |

## Output

Units:
- Length:                 [L] = m
- Mass:                   [m] = g/mol
- Time:                   [t] = s
- Temperature:            [T] = K
- Pressure:               [p] = MPa
- Density:                [ρ] = g/ml
- Energy:                 [E] = eV
- Viscosity:              [η] = Pa*s
- Diffusion Coefficient:  [D] = m²/s
- Thermal Conductivity:   [λ] = W/(m*K)


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

## Input File

__Name__: 'INPUT.txt' (has to be in same folder as 'Master' file)

__Format__:
  - first line: '#' + 'keyword'
  - second line: value assigned to keyword

__Keywords:__

| Name           | Type                        | Description                                                           |
| -------------- | --------------------------- | --------------------------------------------------------------------- |
| modus          | string ["transport","vle"]  | defines the present simulation as 'transport' or 'vle' Simulations    |
| folder         | string                      | path to main folder containing all simulation data of one thermodynamic state (can occur multiple times), asterisk at the end → all folders one level below are evaluated  |
| ensemble       | string ["NVT","NVE","NpT"]  | ensemble to evaluate                                                  |
| timesteps_EQU  | interger (≥ 0)              | number of timesteps to ignore at the start of each simulation         |
| DO_evalulation | integer [0,1]               | 1 - evaluate single folders, 0 - single folders already evaluated     |
| DO_state       | integer [0,1]               | 1 - evaluate complete thermodynamic state (main folder)               |
| N_boot         | integer (≈ No. simulations) | number of bootstrapping repetitions                                   |

Example:

>#folder
>
>F:/MD_Bulk/Others/Methane/2020-04-16_trappe-ua/SIM_T_273.15K_rho_0.1gml
>
>#folder
>
>F:/MD_Bulk/Others/Methane/2020-04-16_trappe-ua/*
>
>#ensemble
>
>NVT
>
>#timesteps_EQU
>
>0
>
>#DO_evaluation
>
>1
>
>#DO_state
>
>1
>
>#N_boot
>
>5

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

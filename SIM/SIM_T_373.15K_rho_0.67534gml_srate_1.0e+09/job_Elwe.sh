#!/bin/bash
# ******* Job name *******
#SBATCH -J SLLOD		 		# name of the job 
#SBATCH -o SIM.%x.%j.out		# Output: %j expands to jobid
#SBATCH -e SIM.%x.%j.err		# Error: %j expands to jobid
# ******* Project to charge *******
#SBATCH --account=TUK-MTD
# ******* Mail notification *******
#SBATCH --mail-type=ALL		
# ******* Run specification ******* 
#SBATCH --nodes=1
#SBATCH --tasks-per-node=24
#SBATCH --time=96:00:00
#SBATCH --mem 10000             
#SBATCH --partition=skylake-96
#SBATCH --licenses=impi

#Empty line above indicates end of SBATCH options
module purge					# clean all module versions
module load intel/2020			# load your version of Intel compilers

mpiexec.hydra -n $SLURM_NTASKS ~/lmp_icc_mpich -in in.Bulk_NVT_sllod

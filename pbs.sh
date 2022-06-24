#!/bin/bash
#SBATCH -A MST110367
#SBATCH -J J_t56_SPE_56_JP
#SBATCH -p dc-MST110367-01              # Partiotion name
#SBATCH --ntasks-per-node=1             # Number of MPI tasks (i.e. processes)
#SBATCH --cpus-per-task=56               # Number of cores per MPI task (OMP)
#SBATCH --nodes=1                       # Maximum number of nodes to be allocated
#SBATCH --mail-type=ALL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jimmy3769tw@gmail.com  # Where to send mail.  Set this to your email address
#SBATCH -t 96:00:00                     # (-t) Wall time limit (days-hrs:min:sec)
#SBATCH --output=%j.log                 # (-o) Path to the standard output and error files relative to the working directory
#SBATCH --error=%j.err                  # (-e) Path to the standard error ouput

module purge
module load compiler/intel/2021   IntelMPI/2021

export OMP_NUM_THREADS=56

export OMP_PLACES=cores

export OMP_PROC_BIND=close
# export OMP_PROC_BIND=spread

export OMP_SCHEDULE=STATIC

# Set OMP_NUM_THREADS to the number of CPUs per task we asked for.
# export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

./mx
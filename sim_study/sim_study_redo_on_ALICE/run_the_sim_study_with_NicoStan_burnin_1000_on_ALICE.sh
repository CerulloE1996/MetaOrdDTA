#!/bin/bash
#SBATCH --job-name=MetaOrdDTA_sim_NicoStan_burnin_1000
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/home/e/ec376/sim_study_runs/logs/%x_%j.out
#SBATCH --partition=short
#SBATCH --qos=normal
## The alpha = (1, ..., 1) sim study of the paper, the cells not yet finished, fitted with NicoStan with 1000 burn-in
## iterations (every other setting the sim's), in the sim's own results folder so that every finished cell is
## skipped. SIM_MODELS picks the models of this job, SIM_N_THREADS the CPUs its workers use (set by sbatch
## --cpus-per-task), SIM_SEED_OFFSET its seeds.
module load R/4.3.1 > /dev/null 2>&1
ulimit -s unlimited 2> /dev/null || ulimit -s 262144
export R_LIBS_USER=$HOME/R_library_4.3.1
export METAORDDTA_OUTPUT_DIR=$HOME/MetaOrdDTA_outputs
export SIM_N_THREADS=$SLURM_CPUS_PER_TASK
export N_BURNIN=1000
export SIM_SEED_OFFSET=${SIM_SEED_OFFSET:-0}
echo "node $(hostname), job $SLURM_JOB_ID, NicoStan burn-in $N_BURNIN, models $SIM_MODELS, CPUs $SIM_N_THREADS, seed offset $SIM_SEED_OFFSET, start $(date)"
Rscript --vanilla $HOME/R_packages/MetaOrdDTA/sim_study/MetaOrdDTA_run_sims_ALICE.R
echo "end $(date), exit $?"

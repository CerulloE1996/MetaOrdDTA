#!/bin/bash
#SBATCH --job-name=whole_sim_MCSE_0_125
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/home/e/ec376/sim_study_runs/logs/%x_%j.out
#SBATCH --partition=short
#SBATCH --qos=normal
## The whole sim study at MCSE 0.125, burn-in 500. SIM_MODELS picks the models of this job, SIM_N_THREADS the CPUs of its workers
## (sbatch --cpus-per-task), SIM_SEED_OFFSET its seeds. Saved seeds are reused, so a job that ends is restarted with the same line:
##   sbatch --cpus-per-task=64 --mem=240G --time=24:00:00 --export=ALL,SIM_MODELS="ord_bivariate" <this file>
module load R/4.3.1 > /dev/null 2>&1
ulimit -s unlimited 2> /dev/null || ulimit -s 262144
export R_LIBS_USER=$HOME/R_library_4.3.1_whole_sim
export METAORDDTA_OUTPUT_DIR=$HOME/MetaOrdDTA_outputs
export SIM_FROZEN_LIBRARY=$HOME/R_library_4.3.1_whole_sim
export SIM_N_THREADS=$SLURM_CPUS_PER_TASK
export SIM_SEED_OFFSET=${SIM_SEED_OFFSET:-0}
export MCSE_TARGET=${MCSE_TARGET:-1}
## (SETTINGS_OF_EVERY_CELL set to "" = the defaults above for every cell)
export SETTINGS_OF_EVERY_CELL=${SETTINGS_OF_EVERY_CELL-$HOME/sim_study_runs/settings_of_every_cell_of_the_whole_sim.csv}
export R_USER_CACHE_DIR=/local/$SLURM_JOB_ID/cache_of_the_job
mkdir -p $R_USER_CACHE_DIR
export N_BURNIN=500
export CLIP_ITER=40
export CLIP_ITER_TAU=120
export LEARNING_RATE=0.05
export TAU_INITIAL=pi
export DIFFUSION_HMC=1
export TAU_JITTER_BURNIN=halton
export SHARE_TAU_II=1
export ESJD_JUMP_POWER=2
export METRIC_START_ITER=85
export METRIC_UPDATE_RESCALES=0
export EPS_REINIT_AT_HANDOVER=1
export INITIAL_EPS_SEARCH_AT_EACH_CHAIN=1
export KAPPA_PRIOR=wide
export MATCHED_PRIORS=band
export SIM_OUTPUT_FOLDER=${SIM_OUTPUT_FOLDER:-whole_sim_MCSE_0.125_NicoStan_burnin_500_alpha_prior_flat_1_band_priors}
echo "node $(hostname), job $SLURM_JOB_ID, models $SIM_MODELS, CPUs $SIM_N_THREADS, seed offset $SIM_SEED_OFFSET, start $(date)"
Rscript --vanilla $HOME/R_packages/MetaOrdDTA/sim_study/MetaOrdDTA_run_sims_ALICE.R
echo "end $(date), exit $?"

#!/bin/bash
#SBATCH --job-name=MetaOrdDTA_sim_NicoStan_band_priors
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/home/e/ec376/sim_study_runs/logs/%x_%j.out
#SBATCH --partition=short
#SBATCH --qos=normal
## The alpha = (1, ..., 1) sim study of the paper, rerun with the priors matched over the whole band of thresholds:
## log(kappa) ~ N(log(150), 1.25) for both random-cutpoint models (KAPPA_PRIOR=wide), the band-matched Jones priors and
## the ord_HSROC log-scale SD of 0.375 for the GAD-2 (MATCHED_PRIORS=band), in its own results folder (SIM_OUTPUT_FOLDER);
## the library with the 0.25 SD inits, burn-in 1000 with clip_iter 70 / clip_iter_tau 200 and learning rate 0.0625
## so that no cell of the old run is reused. SIM_MODELS picks the models of this job, SIM_N_THREADS the CPUs its workers
## use (set by sbatch --cpus-per-task), SIM_SEED_OFFSET its seeds. Submit, e.g.:
##   sbatch --cpus-per-task=64 --mem=240G --time=24:00:00 --export=ALL,SIM_MODELS="ord_bivariate" <this file>
module load R/4.3.1 > /dev/null 2>&1
ulimit -s unlimited 2> /dev/null || ulimit -s 262144
export R_LIBS_USER=$HOME/R_library_4.3.1_newinits
export METAORDDTA_OUTPUT_DIR=$HOME/MetaOrdDTA_outputs
export SIM_N_THREADS=$SLURM_CPUS_PER_TASK
export N_BURNIN=1000
export CLIP_ITER=70
export CLIP_ITER_TAU=200
export LEARNING_RATE=0.0625
export SIM_SEED_OFFSET=${SIM_SEED_OFFSET:-0}
export KAPPA_PRIOR=wide
export MATCHED_PRIORS=band
export SIM_OUTPUT_FOLDER=${SIM_OUTPUT_FOLDER:-revised_DGMs_fixed_raw_C_J_sim_results_NicoStan_MCSE_0.125_alpha_prior_flat_1_band_priors}
echo "node $(hostname), job $SLURM_JOB_ID, NicoStan burn-in $N_BURNIN, models $SIM_MODELS, CPUs $SIM_N_THREADS, seed offset $SIM_SEED_OFFSET, clip $CLIP_ITER / $CLIP_ITER_TAU, LR $LEARNING_RATE, library $R_LIBS_USER, kappa prior $KAPPA_PRIOR, matched priors $MATCHED_PRIORS, folder $SIM_OUTPUT_FOLDER, start $(date)"
Rscript --vanilla $HOME/R_packages/MetaOrdDTA/sim_study/MetaOrdDTA_run_sims_ALICE.R
echo "end $(date), exit $?"

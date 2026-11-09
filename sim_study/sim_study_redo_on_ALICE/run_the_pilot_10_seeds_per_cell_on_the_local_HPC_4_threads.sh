#!/bin/bash
## The 10-seeds-per-cell pilot of the sim study on the local HPC, 4 threads on CPUs 180-191. Every cell runs the same
## sampler settings (burn-in 500, clip_iter 70, clip_iter_tau 200, learning rate 0.0625, every other setting the
## runner's own default) and 1000 iterations; no per-cell settings. The ordinal-bivariate (fixed and random cutpoints)
## and R&G models first, the stratified bivariate last. Start it with:
##   systemd-run --user --unit=pilot_10_seeds_HPC_4_threads -p MemoryMax=40G -p MemorySwapMax=0 \
##     bash -c "taskset -c <4 CPUs> bash <this file> > <log> 2>&1"
ulimit -s unlimited 2> /dev/null || ulimit -s 262144
export SIM_PACKAGE_DIR=/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA
export SIM_FROZEN_LIBRARY=/home/enzocerullo/Documents/Work/PhD_work/R_packages/R_library_for_the_sim_pilot_on_the_local_HPC
export R_LIBS_USER=$SIM_FROZEN_LIBRARY
export METAORDDTA_OUTPUT_DIR=/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA_outputs
export SIM_N_THREADS=${SIM_N_THREADS:-4}
export SIM_N_SEEDS=${SIM_N_SEEDS:-10}
export SIM_SEED_OFFSET=${SIM_SEED_OFFSET:-0}
export MCSE_TARGET=0.125
export SIM_DGMS=${SIM_DGMS:-ord_biv_random_C,R&G_random_C,ord_biv_fixed_C}
export SIM_OUTPUT_FOLDER=${SIM_OUTPUT_FOLDER:-pilot_10_seeds_per_cell_one_setting_for_every_cell}
export SAMPLING_ALGORITHM=${SAMPLING_ALGORITHM:-NicoStan}
## (no per-cell settings: every cell runs the settings below)
export SETTINGS_OF_EVERY_CELL=""
export N_ITER_FOR_4_CHAINS=1000
export N_BURNIN=500
export CLIP_ITER=70
export CLIP_ITER_TAU=200
export LEARNING_RATE=0.0625
driver=$SIM_PACKAGE_DIR/sim_study/sim_study_redo_on_ALICE/MetaOrdDTA_run_sims_ALICE.R
echo "local HPC, sampler $SAMPLING_ALGORITHM, DGMs $SIM_DGMS, CPUs $SIM_N_THREADS, seeds $SIM_N_SEEDS from $SIM_SEED_OFFSET, start $(date)"
SIM_MODELS=ord_bivariate SIM_RANDOM_THRESHOLDS=FALSE,TRUE Rscript --vanilla $driver
SIM_MODELS="R&G" SIM_RANDOM_THRESHOLDS=TRUE Rscript --vanilla $driver
SIM_MODELS=Bivariate SIM_RANDOM_THRESHOLDS=FALSE Rscript --vanilla $driver
echo "end $(date), exit $?"

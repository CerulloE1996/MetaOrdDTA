#!/bin/bash
## The 10-seeds-per-cell pilot of the sim study on the local HPC: 8 threads on CPUs 180-191, the sampler settings of
## run_the_whole_sim_at_MCSE_0_125_on_ALICE.sh. SIM_DGMS picks the data-generating models of this run, SIM_MODELS its
## models. Start it with:
##   systemd-run --user --unit=pilot_10_seeds_HPC_8_threads -p MemoryMax=40G -p MemorySwapMax=0 \
##     bash -c "taskset -c 180,183,184,185,186,187,188,190 bash <this file> > <log> 2>&1"
ulimit -s unlimited 2> /dev/null || ulimit -s 262144
export SIM_PACKAGE_DIR=/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA
export SIM_FROZEN_LIBRARY=/home/enzocerullo/Documents/Work/PhD_work/R_packages/R_library_for_the_sim_pilot_on_the_local_HPC
export R_LIBS_USER=$SIM_FROZEN_LIBRARY
export METAORDDTA_OUTPUT_DIR=/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA_outputs
export SIM_N_THREADS=${SIM_N_THREADS:-8}
export SIM_N_SEEDS=${SIM_N_SEEDS:-10}
export SIM_SEED_OFFSET=${SIM_SEED_OFFSET:-0}
export MCSE_TARGET=0.125
export SIM_MODELS=${SIM_MODELS:-ord_bivariate}
export SIM_DGMS=${SIM_DGMS:-Jones}
export SIM_OUTPUT_FOLDER=${SIM_OUTPUT_FOLDER:-pilot_10_seeds_per_cell_prior_mean_anchor_9_Nov}
export SAMPLING_ALGORITHM=${SAMPLING_ALGORITHM:-NicoStan}
## (SETTINGS_OF_EVERY_CELL set to "" = the defaults below for every cell)
export SETTINGS_OF_EVERY_CELL=${SETTINGS_OF_EVERY_CELL-$SIM_PACKAGE_DIR/sim_study/sim_study_redo_on_ALICE/settings_of_every_cell_of_the_whole_sim.csv}
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
echo "local HPC, sampler $SAMPLING_ALGORITHM, models $SIM_MODELS, DGMs $SIM_DGMS, CPUs $SIM_N_THREADS, seeds $SIM_N_SEEDS from $SIM_SEED_OFFSET, start $(date)"
Rscript --vanilla $SIM_PACKAGE_DIR/sim_study/sim_study_redo_on_ALICE/MetaOrdDTA_run_sims_ALICE.R
echo "end $(date), exit $?"

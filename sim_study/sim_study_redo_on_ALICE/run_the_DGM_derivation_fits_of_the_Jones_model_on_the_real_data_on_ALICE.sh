#!/bin/bash
#SBATCH --job-name=DGM_derivation_Jones
#SBATCH --partition=short
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=06:00:00
#SBATCH --output=/home/e/ec376/sim_study_runs/logs/%x_%j.out
##
## The DGM derivation, sixth pass, on the cleaned real data (per-test files cut from
## Klaus_collab_NMA_work/derived_data/Klaus_et_al_data/real_data_8_tests_GAD: test 1 = GAD-2, 3 = HADS,
## 4 = BAI), one test per job (sbatch --export=ALL,SIM_TEST_TEST=<GAD_2|HADS|BAI> this file): the sim runner
## fits the Jones model with its converted Se / Sp prior of set 0 (the trivariate Normal of location, log scale and lambda) to the real
## data with NicoStan (the library R_library_4.3.1_Jones_converted), 1000 burn-in (clip 70 / 200, learning rate 0.0625,
## the 0.25 SD inits), 4000 sampling iterations on each of 16 chains; one test per job
## (sbatch --export=ALL,SIM_TEST_TEST=<GAD_2|HADS|BAI> this file). Each job writes its own TEST_
## folder, so the three jobs run at once.
module load R/4.3.1 > /dev/null 2>&1
ulimit -s unlimited 2> /dev/null || ulimit -s 262144
export R_LIBS_USER=$HOME/R_library_4.3.1_Jones_converted
export METAORDDTA_OUTPUT_DIR=$HOME/MetaOrdDTA_outputs
export SIM_TEST=1
export SIM_N_SEEDS=1
export SIM_SEED_OFFSET=0
export SIM_N_THREADS=16
export SIM_N_CHAINS_PER_WORKER=16
export SIM_MODELS="Jones"
export SIM_RUN_LABEL=${SIM_RUN_LABEL:-_Jones}
export SIM_REAL_DATA_DIR=$HOME/sim_study_runs/real_data_for_the_DGM_derivation_cleaned
export SIM_TEST_N_STUDIES=50
export SIM_TEST_DGM="R&G_random_C"
export N_ITER_FOR_4_CHAINS=16000
export N_BURNIN=1000
export BURNIN_ALGORITHM=ESJD
export CLIP_ITER=70
export CLIP_ITER_TAU=200
export LEARNING_RATE=0.0625
export KAPPA_PRIOR=wide
export MATCHED_PRIORS=band
export SIM_OUTPUT_FOLDER=DGM_derivation_v6_${SIM_TEST_TEST}${SIM_RUN_LABEL}
T=$HOME/MetaOrdDTA_outputs/sim_study/TEST_${SIM_OUTPUT_FOLDER}
OUT=$HOME/MetaOrdDTA_outputs/sim_study/DGM_derivation_v6
mkdir -p $OUT
[ -d $T ] && mv $T ${T}_set_aside_$SLURM_JOB_ID
echo "== test $SIM_TEST_TEST on NicoStan, job $SLURM_JOB_ID, start $(date)"
Rscript --vanilla $HOME/R_packages/MetaOrdDTA/sim_study/MetaOrdDTA_run_sims_ALICE.R \
      > $HOME/sim_study_runs/logs/DGM_derivation_v6_${SIM_TEST_TEST}_$SLURM_JOB_ID.log 2>&1
[ -d $T ] && mv $T $OUT/${SIM_TEST_TEST}${SIM_RUN_LABEL}_NicoStan
for g in $OUT/${SIM_TEST_TEST}${SIM_RUN_LABEL}_NicoStan/gate_of_every_fit*; do echo "$g" | sed "s|.*gate_of_every_fit_||" | cut -c1-110; tail -n +2 "$g" | cut -d, -f1-7; done
grep -h "ERROR" $OUT/${SIM_TEST_TEST}${SIM_RUN_LABEL}_NicoStan/log_of_the_workers.txt 2> /dev/null | cut -c1-200
echo "== test $SIM_TEST_TEST end $(date)"

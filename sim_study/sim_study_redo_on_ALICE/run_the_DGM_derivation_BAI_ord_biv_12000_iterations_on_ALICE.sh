#!/bin/bash
#SBATCH --job-name=DGM_derivation_v5_BAI_ordbiv_12000
#SBATCH --partition=short
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=120G
#SBATCH --time=06:00:00
#SBATCH --output=/home/e/ec376/sim_study_runs/logs/%x_%j.out
##
## The DGM re-derivation, fifth pass, on the CLEANED real data (per-test files cut from
## Klaus_collab_NMA_work/derived_data/Klaus_et_al_data/real_data_8_tests_GAD: test 1 = GAD-2, 3 = HADS, 4 = BAI), one
## (the BAI ord-biv fits again with 12000 sampling iterations over 4 chains: the random-cutpoint fit of the 4000-iteration run had a tail ESS of 633 and an SD ESS of 952 on one Se / Sp, under the 1000 rule)
## fixed-C, ord-biv random-C pop-median, R&G random-C pop-median; the final flat alpha = 1 priors, kappa prior C, the
## band-matched priors) to the real data with NicoStan, 1000 burn-in (clip 70 / 200, learning rate 0.0625, the 0.25 SD
## inits), 4000 sampling iterations over 4 chains. Each job writes its own TEST_ folder, so the three jobs run at once.
module load R/4.3.1 > /dev/null 2>&1
ulimit -s unlimited 2> /dev/null || ulimit -s 262144
export R_LIBS_USER=$HOME/R_library_4.3.1_newinits:$HOME/R_library_4.3.1
export METAORDDTA_OUTPUT_DIR=$HOME/MetaOrdDTA_outputs
export SIM_TEST=1
export SIM_N_SEEDS=1
export SIM_SEED_OFFSET=0
export SIM_N_THREADS=12
export SIM_MODELS="ord_bivariate"
export SIM_REAL_DATA_DIR=$HOME/sim_study_runs/real_data_for_the_DGM_derivation_cleaned
export SIM_TEST_N_STUDIES=50
export SIM_TEST_DGM="R&G_random_C"
export N_ITER_FOR_4_CHAINS=12000
export N_BURNIN=1000
export BURNIN_ALGORITHM=ESJD
export CLIP_ITER=70
export CLIP_ITER_TAU=200
export LEARNING_RATE=0.0625
export KAPPA_PRIOR=wide
export MATCHED_PRIORS=band
export SIM_OUTPUT_FOLDER=DGM_derivation_v5_${SIM_TEST_TEST}_ordbiv_12000
T=$HOME/MetaOrdDTA_outputs/sim_study/TEST_${SIM_OUTPUT_FOLDER}
OUT=$HOME/MetaOrdDTA_outputs/sim_study/DGM_derivation_v5
mkdir -p $OUT
[ -d $T ] && mv $T ${T}_set_aside_$SLURM_JOB_ID
echo "== test $SIM_TEST_TEST on NicoStan, job $SLURM_JOB_ID, start $(date)"
Rscript --vanilla $HOME/R_packages/MetaOrdDTA/sim_study/MetaOrdDTA_run_sims_ALICE.R \
      > $HOME/sim_study_runs/logs/DGM_derivation_v5_${SIM_TEST_TEST}_$SLURM_JOB_ID.log 2>&1
[ -d $T ] && mv $T $OUT/${SIM_TEST_TEST}_NicoStan_ordbiv_12000
for g in $OUT/${SIM_TEST_TEST}_NicoStan_ordbiv_12000/gate_of_every_fit*; do echo "$g" | sed "s|.*gate_of_every_fit_||" | cut -c1-110; tail -n +2 "$g" | cut -d, -f1-7; done
grep -h "ERROR" $OUT/${SIM_TEST_TEST}_NicoStan_ordbiv_12000/log_of_the_workers.txt 2> /dev/null | cut -c1-200
echo "== test $SIM_TEST_TEST end $(date)"

#!/bin/bash
#SBATCH --job-name=band_prior_fits
#SBATCH --partition=short
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=08:00:00
#SBATCH --output=/home/e/ec376/sim_study_runs/matched_priors_band/logs/%x_%A_%a.out
##
## The NicoStan prior-only fits, matched on the whole band of thresholds of each test, of the list file (one line
## per fit: model test SD1 SD2 SD3 kappa seed [mean1 mean2 [mean3]]),
## this array task taking every N_TASKS-th line from its index, 16 fits at once on the node (4 chains per fit). With
## COMPILE_ONLY=1 it runs the first fit of each model alone, so that every model is built once before the array
## starts. LIST names the list file (default list_of_the_fits.txt).
module load R/4.3.1 > /dev/null 2>&1
ulimit -s unlimited 2> /dev/null || ulimit -s 262144
export R_LIBS_USER=${R_LIBS_OF_THE_FITS:-$HOME/R_library_4.3.1}
export METAORDDTA_OUTPUT_DIR=$HOME/MetaOrdDTA_outputs
D=$HOME/R_packages/MetaOrdDTA/sim_study
ROOT=$HOME/sim_study_runs/matched_priors_band
LIST=${LIST:-$ROOT/list_of_the_band_fits.txt}
N_CHAINS=4; N_WARMUP=250; N_ITER=500
export D ROOT N_CHAINS N_WARMUP N_ITER
run_one_fit() {
      model=$1; test=$2; a=$3; b=$4; c=$5; kappa=$6; seed=$7; m1=${8:-}; m2=${9:-}; m3=${10:-}
      OUT=$ROOT/${model}_${test}_kappa_${kappa}
      mkdir -p $OUT/logs
      MEANS=""; MEAN_ARGS=""
      if [ -n "$m1" ]; then MEANS="_means_${m1}_${m2}"; MEAN_ARGS="$m1 $m2"; fi
      if [ -n "$m3" ]; then MEANS="${MEANS}_${m3}"; MEAN_ARGS="$MEAN_ARGS $m3"; fi
      KAPPA_PART="_kappa_${kappa}"
      case $model in jones|xu_fixed) KAPPA_PART="";; esac
      FILE=$OUT/${model}_${test}_SDs_${a}_${b}_${c}${MEANS}${KAPPA_PART}_seed_${seed}.rds
      [ -f $FILE ] && return 0
      KAPPA_PRIOR=$kappa Rscript $D/run_one_NicoStan_PO_fit_alpha_1_matched_priors.R \
            $model $test $N_CHAINS $N_WARMUP $N_ITER $OUT $seed $a $b $c $MEAN_ARGS \
            > $OUT/logs/SDs_${a}_${b}_${c}${MEANS}_seed_${seed}.log 2>&1
}
export -f run_one_fit
START=$(date +%s)
echo "node $(hostname), job $SLURM_JOB_ID task ${SLURM_ARRAY_TASK_ID:-compile}, list $LIST, start $(date)"
if [ "${COMPILE_ONLY:-0}" = "1" ]; then
      for m in xu_fixed xu_rand_popmed rg_rand_popmed jones; do
            line=$(awk -v m=$m '$1 == m' $LIST | head -1)
            [ -z "$line" ] && continue
            run_one_fit $line
            echo "$m built and its first fit done after $(( $(date +%s) - START )) s: $(grep -h "^DONE\|Error" $ROOT/${m}_*/logs/*.log | head -1 | cut -c1-200)"
      done
      exit 0
fi
## (the first fit of the task alone, so that the node's NicoStan runtime and the model builds are made once)
FIRST=$(awk -v i=$SLURM_ARRAY_TASK_ID -v n=$N_TASKS 'NR % n == i' $LIST | head -1)
run_one_fit $FIRST
awk -v i=$SLURM_ARRAY_TASK_ID -v n=$N_TASKS 'NR % n == i' $LIST | tail -n +2 | \
      xargs -P 16 -L 1 bash -c 'run_one_fit $0 $1 $2 $3 $4 $5 $6 $7 $8 $9'
echo "task $SLURM_ARRAY_TASK_ID done after $(( $(date +%s) - START )) s"

#!/bin/bash
#SBATCH --job-name=four_hardest_cells_500_40_120_newinits
#SBATCH --partition=short
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --time=04:00:00
#SBATCH --output=/home/e/ec376/sim_study_runs/logs/%x_%j.out
##
## The ord_HSROC random-cutpoint pop-median model on the BAI R&G (10 and 50 studies) and HADS R&G data, the ord_biv on BAI and
## The 4 hardest cells, 30 seeds each, 500 burn-in iterations with clip_iter 40 / clip_iter_tau 120 and learning rate
## 0.0625, on the library with the 0.25 SD inits; the burn-in record of every fit saved.
module load R/4.3.1 > /dev/null 2>&1
ulimit -s unlimited 2> /dev/null || ulimit -s 262144
export R_LIBS_USER=$HOME/R_library_4.3.1_newinits
export METAORDDTA_OUTPUT_DIR=$HOME/MetaOrdDTA_outputs
export SIM_TEST=1
export SIM_SEED_OFFSET=0
export SIM_N_SEEDS=30
export SIM_N_THREADS=64
export DUMP_BURNIN_TRACE=1
export SIM_MODELS="R&G"
T=$HOME/MetaOrdDTA_outputs/sim_study/TEST_fixed_raw_C_J_sim_results_NicoStan_MCSE_0.125_alpha_prior_flat_1
[ -d $T ] && mv $T ${T}_set_aside_before_the_500_test
run_variant() {
      name=$1
      shift
      echo "== setting $name ($*), start $(date)"
      env "$@" Rscript --vanilla $HOME/R_packages/MetaOrdDTA/sim_study/MetaOrdDTA_run_sims_ALICE.R \
            > $HOME/sim_study_runs/logs/ord_HSROC_arm_${name}_$SLURM_JOB_ID.log 2>&1
      [ -d $T ] && mv $T ${T}_RG_30_seeds_newinits_${name}
      g=$(ls ${T}_RG_30_seeds_newinits_${name}/gate_of_every_fit* 2> /dev/null | head -1)
      h=$(ls ${T}_RG_30_seeds_newinits_${name}/HMC_of_every_fit* 2> /dev/null | head -1)
      if [ -n "$g" ]; then
            tail -n +2 "$g" | awk -F, '{n++; k += ($7 == "TRUE"); d += $6; if ($6 > 1) {m++; dm += $6}} END {printf "   %d fits, %d kept, %d over 1%% divergent (mean %.1f%% among them), mean divergent %.2f%%\n", n, k, m, (m > 0 ? dm / m : 0), d / n}'
            echo "   seed | min_ESS | pct_div | kept | s_burnin | s_sampling | n_grads | ESS_per_1000_grads | eps | tau | L_mean"
            join -t, -j1 <(tail -n +2 "$g" | cut -d, -f1,2,6,7,8,9,10,11 | sort -t, -k1,1) <(tail -n +2 "$h" | cut -d, -f1,4,5,6 | sort -t, -k1,1) | sort -t, -k1,1n | awk -F, '{printf "   %s | %.0f | %.1f | %s | %.0f | %.1f | %d | %.2f | %.2e | %.1f | %.0f\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}'
      else
            echo "   no gate file: tail of the log"
            tail -n 25 $HOME/sim_study_runs/logs/ord_HSROC_arm_${name}_$SLURM_JOB_ID.log | cut -c1-160
      fi
      echo "== setting $name end $(date)"
}
run_variant b500_BAI_RGdgm_n10_clip_40_tau_120_LR_0_0625 SIM_TEST_TEST=BAI SIM_TEST_DGM=R\&G_random_C SIM_TEST_N_STUDIES=10 BURNIN_ALGORITHM=ESJD N_BURNIN=500 CLIP_ITER=40 CLIP_ITER_TAU=120 LEARNING_RATE=0.0625
run_variant b500_HADS_RGdgm_n10_clip_40_tau_120_LR_0_0625 SIM_TEST_TEST=HADS SIM_TEST_DGM=R\&G_random_C SIM_TEST_N_STUDIES=10 BURNIN_ALGORITHM=ESJD N_BURNIN=500 CLIP_ITER=40 CLIP_ITER_TAU=120 LEARNING_RATE=0.0625
run_variant b500_BAI_OBRCdgm_n10_clip_40_tau_120_LR_0_0625 SIM_TEST_TEST=BAI SIM_TEST_DGM=ord_biv_random_C SIM_TEST_N_STUDIES=10 BURNIN_ALGORITHM=ESJD N_BURNIN=500 CLIP_ITER=40 CLIP_ITER_TAU=120 LEARNING_RATE=0.0625
run_variant b500_BAI_RGdgm_n50_clip_40_tau_120_LR_0_0625 SIM_TEST_TEST=BAI SIM_TEST_DGM=R\&G_random_C SIM_TEST_N_STUDIES=50 BURNIN_ALGORITHM=ESJD N_BURNIN=500 CLIP_ITER=40 CLIP_ITER_TAU=120 LEARNING_RATE=0.0625
run_variant b500_HADS_OBRCdgm_n10_s21_30_clip_40_tau_120_LR_0_0625 SIM_SEED_OFFSET=20 SIM_N_SEEDS=10 SIM_TEST_TEST=HADS SIM_TEST_DGM=ord_biv_random_C SIM_TEST_N_STUDIES=10 BURNIN_ALGORITHM=ESJD N_BURNIN=500 CLIP_ITER=40 CLIP_ITER_TAU=120 LEARNING_RATE=0.0625

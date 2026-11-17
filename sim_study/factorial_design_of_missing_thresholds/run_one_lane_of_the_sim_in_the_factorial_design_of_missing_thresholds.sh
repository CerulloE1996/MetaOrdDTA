#!/bin/bash
## One lane of the sim of the MetaOrdDTA models in the factorial design of missing thresholds (every test at 15, 55 and 85% of its
## thresholds missing): the cells of a list one after another on 4 CPUs, through the sim's driver and runner, each cell until the
## MCSE target of 0.125 (counted over the fits of every lane) or until the number of seeds given is fitted; a seed already fitted is
## never fitted again. Fused nico::math file, NicoStan, the sim's sampler settings, 4 chains, 500 burn-in; the draws are not kept.
## Every setting of a cell is written into ONE settings file (kept in the folder of the logs) that the driver reads; the driver
## and the runner read no environment variable. The package folder and the library of a cell are given with the cell; the lane's
## cache holds the compiled model of every library. Results (TEST_<output folder>_with_<percent>_pct_of_thresholds_missing, one
## folder for each proportion missing) and the cache are on the home drive.
## Usage: bash <this file> <list file> <first of its 4 CPUs> <number of this lane, for its ports> <seed offset> <number of seeds a cell at most>
##                         <output folder, without TEST_> <folder of the logs> <cache folder of this lane>
##                         [<.rds file of sampler settings that replace the runner's own; not given = the runner's own>]
##   The list file, one cell per line: <model: ord_bivariate | ord_HSROC | Bivariate | Jones> <random cutpoints: TRUE | FALSE> <test>
##                                     <number of studies> <DGM> <prior set> <sampling iterations per chain> <package folder> <library>
##                                     <proportion of the thresholds missing: 0.15 | 0.55 | 0.85>
[ $# -eq 8 ] || [ $# -eq 9 ] || { echo "eight or nine arguments are needed"; exit 1; }
list_file=$1
first_cpu=$2
number=$3
seed_offset=$4
n_seeds=$5
output_folder=$6
logs=$7
cache=$8
sampler_settings_file=NULL
[ $# -eq 9 ] && { [ -f "$9" ] || { echo "the file of sampler settings is not there: $9"; exit 1; }; sampler_settings_file="\"$9\""; }
W=/home/enzocerullo/Documents/Work/PhD_work/R_packages
mkdir -p $logs
case "$cache" in /home/enzocerullo/.cache/*) ;; *) echo "the cache of a lane must be on the home drive"; exit 1;; esac
[ -d "$cache/R" ] || { echo "the cache of the lane has no compiled models"; exit 1; }
export PATH=/home/enzocerullo/R/R-4.3.3/bin:$PATH
export NUISANCE_TO_TRACK=none
export R_USER_CACHE_DIR=$cache
while read -r model random_cutpoints test n_studies DGM prior_set n_iter P L proportion_missing; do
      [ -z "$model" ] && continue
      [ -f "$P/sim_study/run_MULTIPLE_MetaOrd_simulation.R" ] || { echo "the package folder has no runner: $P"; exit 1; }
      [ -d "$L/MetaOrdDTA" ] || { echo "the library has no MetaOrdDTA: $L"; exit 1; }
      case "$proportion_missing" in 0.15|0.55|0.85) ;; *) echo "the proportion missing of a cell must be 0.15, 0.55 or 0.85: $proportion_missing"; exit 1;; esac
      percent_missing=${proportion_missing#0.}
      export R_LIBS_USER=$L
      driver=$P/sim_study/sim_study_redo_on_ALICE/MetaOrdDTA_run_sims_ALICE.R
      name="lane_${number}_${model}_random_cutpoints_${random_cutpoints}_${test}_${n_studies}_studies_DGM_${DGM}_${prior_set}_${percent_missing}_pct_missing"
      log="$logs/log_of_$name.txt"
      settings="$logs/settings_of_$name.R"
      {
            echo "## The settings of one run of the sim driver: every setting of the run, each written once."
            echo "frozen_library <- \"$L\""
            echo "package_dir <- \"$P\""
            echo "results_root <- \"$W/MetaOrdDTA_outputs\""
            echo "output_folder <- \"${output_folder}_with_${percent_missing}_pct_of_thresholds_missing\""
            echo "results_in_a_TEST_folder <- TRUE"
            echo "n_threads <- 4"
            echo "n_chains_per_worker <- 4"
            echo "port_start <- $((14000 + 40 * number))"
            echo "seed_offset <- $seed_offset"
            echo "n_seeds_per_cell <- $n_seeds"
            echo "MCSE_target <- 0.125"
            echo "sampling_algorithm <- \"NicoStan\""
            echo "DGMs <- \"$DGM\""
            echo "models <- \"$model\""
            echo "random_thresholds <- $random_cutpoints"
            echo "prior_sets <- \"$prior_set\""
            echo "tests <- \"$test\""
            echo "n_studies <- $n_studies"
            echo "pop_median_models <- TRUE"
            echo "Stan_math_only_model <- FALSE"
            echo "real_data_dir <- NULL"
            echo "n_burnin <- 500"
            echo "n_iter_for_4_chains <- $n_iter"
            echo "log_lik_indicator <- 0"
            echo "nuisance_params <- NULL"
            echo "keep_the_data_file_in <- NULL"
            echo "keep_the_divergence_trace_in <- NULL"
            echo "sampler_settings_to_test_file <- $sampler_settings_file"
            echo "target_proportion_missing <- $proportion_missing"
      } > "$settings"
      echo "$(date '+%d %b %H:%M:%S') start: $model (random cutpoints $random_cutpoints), $test, $n_studies studies, DGM $DGM, $prior_set, $percent_missing% missing, $n_iter iterations, CPUs $first_cpu-$((first_cpu + 3))"
      echo "kind fused, model $model, test $test, studies $n_studies, DGM $DGM, prior set $prior_set, thresholds missing $percent_missing%, seeds $((seed_offset + 1)) to $((seed_offset + n_seeds)), CPUs $first_cpu-$((first_cpu + 3)), random cutpoints $random_cutpoints, settings file $settings, library $L, package folder $P, start $(date)" >> "$log"
      taskset -c $first_cpu-$((first_cpu + 3)) Rscript --vanilla $driver "$settings" >> "$log" 2>&1 < /dev/null
      echo "$(date '+%d %b %H:%M:%S') end:   $model (random cutpoints $random_cutpoints), $test, $n_studies studies, DGM $DGM, $prior_set, $percent_missing% missing, exit $?"
done < "$list_file"
echo "$(date '+%d %b %H:%M:%S') every cell of this list is done"

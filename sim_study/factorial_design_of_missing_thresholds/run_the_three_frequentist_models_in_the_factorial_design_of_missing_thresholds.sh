#!/bin/bash
## The three frequentist models of the sim (Hoyer 2018 time-to-event and the logit GLMM together; the diagmeta model) on the cells
## of the factorial design of missing thresholds: for each proportion missing (55, 85, 15%) the cells of its list one after
## another, each until the MCSE target of 0.125 or the largest number of seeds. Each seed's data set is the one the sim's lanes
## fit for that cell and seed. Results: one folder for each kind of model and proportion missing.
## Usage: bash <this file> <competing | diagmeta> <folder of the lists> <CPUs, e.g. 180-185> <number of cores> <largest number of seeds of a cell>
##                         <folder of the logs>
[ $# -eq 6 ] || { echo "six arguments are needed"; exit 1; }
kind=$1; lists=$2; cpus=$3; n_cores=$4; n_seeds=$5; logs=$6
O=/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA_outputs/sim_study
export PATH=/home/enzocerullo/R/R-4.3.3/bin:$PATH
mkdir -p $logs
for percent in 55 85 15; do
      list=$lists/list_of_the_cells_with_${percent}_pct_missing.txt
      [ -f $list ] || { echo "no list: $list"; exit 1; }
      if [ $kind = competing ]; then
            folder=$O/sim_of_the_competing_models_Hoyer_2018_and_logit_GLMM_in_the_factorial_design_15_Nov_with_${percent}_pct_of_thresholds_missing
            mkdir -p $folder
            echo "$(date '+%d %b %H:%M:%S') start: Hoyer 2018 and logit GLMM, $percent% missing"
            taskset -c $cpus Rscript --vanilla $O/fit_the_competing_models_to_the_cells_of_a_list_until_the_MCSE_target.R $list $n_cores $folder $n_seeds \
                  >> $logs/log_of_Hoyer_2018_and_logit_GLMM_with_${percent}_pct_missing.txt 2>&1 < /dev/null
      elif [ $kind = diagmeta ]; then
            folder=$O/sim_of_the_diagmeta_model_DIDS_logistic_in_the_factorial_design_15_Nov_with_${percent}_pct_of_thresholds_missing
            mkdir -p $folder
            echo "$(date '+%d %b %H:%M:%S') start: diagmeta (DIDS, logistic), $percent% missing"
            taskset -c $cpus Rscript --vanilla $O/fit_the_diagmeta_model_to_the_cells_of_a_list_until_the_MCSE_target.R $list $n_cores $folder $n_seeds 200 DIDS logistic \
                  >> $logs/log_of_the_diagmeta_model_with_${percent}_pct_missing.txt 2>&1 < /dev/null
      else
            echo "the kind must be competing or diagmeta"; exit 1
      fi
      echo "$(date '+%d %b %H:%M:%S') end:   $kind, $percent% missing, exit $?"
done
echo "$(date '+%d %b %H:%M:%S') every proportion missing is done"

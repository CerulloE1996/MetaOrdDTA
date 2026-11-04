#!/bin/bash
#SBATCH --job-name=score_band_new
#SBATCH --partition=short
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/home/e/ec376/sim_study_runs/matched_priors_band/logs/%x_%j.out
module load R/4.3.1 > /dev/null 2>&1
export R_LIBS_USER=$HOME/R_library_4.3.1
Rscript $HOME/R_packages/MetaOrdDTA/sim_study/score_the_matched_prior_fits_on_the_band_against_the_ord_biv_fixed_prior.R \
      $HOME/sim_study_runs/matched_priors_band 0.05

#!/bin/bash
#SBATCH --job-name=install_MetaOrdDTA_Jones_converted
#SBATCH --partition=devel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --chdir=/home/e/ec376
#SBATCH --output=/home/e/ec376/sim_study_runs/logs/install_MetaOrdDTA_Jones_converted_%j.log
#SBATCH --export=NONE
## Installs MetaOrdDTA from $HOME/R_packages/MetaOrdDTA (the synced source) into R_library_4.3.1_Jones_converted, beside
## the NicoStan and NicoSIMD of R_library_4.3.1_whole_sim copied there once: the library of the Jones real-data fits and
## the Jones pilot with the converted Se / Sp prior.
module load R/4.3.1
L=$HOME/R_library_4.3.1_Jones_converted
mkdir -p $L
for p in NicoStan NicoSIMD; do [ -d $L/$p ] || cp -r $HOME/R_library_4.3.1_whole_sim/$p $L/; done
export R_LIBS=$L
export R_LIBS_USER=$L
export MAKEFLAGS=-j4
R CMD INSTALL --no-test-load -l $L $HOME/R_packages/MetaOrdDTA
echo "exit of the install: $?"
grep -c "stan_external_hyperbolic_functions.hpp" $L/NicoSIMD/include/NicoSIMD/stan_external_functions_of_NicoSIMD.hpp

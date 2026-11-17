#!/bin/bash
## After the 17 cells of prior set 0 that had no score are fitted: (1) every O-biv-RC and stratified bivariate cell of prior set 0
## of HADS at 15% and 55% missing and of BAI at 85% missing is scored again with the sim's own script (same rule: a threshold counts
## when at least 1 study reports it; target 0.125; 30 kept fits at least), three small runs at the lowest priority; (2) their rows
## replace the rows of the same cells in the table of the factorial sim; (3) the table of the nine sets, the results tables, the
## ranking tables and the figures are made again. Usage: bash <this file>
S=/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA_outputs/sim_study
U2=$S/update_of_the_paper_with_the_factorial_sim
F=$S/final_results_of_the_sim_in_the_factorial_design_scored_as_it_stands_17_Nov
R=/home/enzocerullo/R/R-4.3.3/bin/Rscript
W=$U2/rescoring_of_the_cells_fitted_after_the_stop
mkdir -p $W $U2/logs
models='param_(ord_bivariatesoft_FALSErand_thr_TRUE|Bivariatesoft_FALSErand_thr_FALSE)set_0'
for p in 15 55 85; do
      test=HADS; [ $p = 85 ] && test=BAI
      printf '%s\t%s\n' "$S/TEST_sim_in_the_factorial_design_of_missing_thresholds_15_Nov_with_${p}_pct_of_thresholds_missing" "test_${test}red_fact.*${models}" > $W/list_of_the_folder_${p}_pct_missing.txt
      ( cd $W && nice -n 19 ionice -c 3 $R --vanilla $S/RMSE_bias_coverage_and_MCSE_of_every_cell_of_every_model_of_the_whole_sim_from_a_list_of_folders.R \
            /home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA $W/list_of_the_folder_${p}_pct_missing.txt $W/table_of_the_cells_scored_again_${p}_pct_missing.rds 4 1 0.125 30 \
            > $W/log_of_the_scoring_${p}_pct_missing.txt 2>&1 ) &
done
wait
grep -h 'results files listed' $W/log_of_the_scoring_*_pct_missing.txt | cut -c1-220
$R --vanilla -e "
F <- '$F'; W <- '$W'
t <- readRDS(file.path(F, 'table_of_every_cell_all_three_shares_missing.rds'))
if (!file.exists(file.path(F, 'table_of_every_cell_all_three_shares_missing_as_scored_at_the_stop.rds'))) saveRDS(t, file.path(F, 'table_of_every_cell_all_three_shares_missing_as_scored_at_the_stop.rds'))
key <- function(d) paste(d\$model, d\$prior_set, d\$test, d\$pct, d\$n_studies, d\$DGM)
for (p in c(15, 55, 85)) {
  n <- readRDS(file.path(W, sprintf('table_of_the_cells_scored_again_%d_pct_missing.rds', p))); n\$pct <- p
  stopifnot(setequal(names(n), names(t)))
  cat(p, '% missing: cells scored again', nrow(n), '; new cells', sum(!key(n) %in% key(t)), '; at the target', sum(n\$mcse_met), '\n')
  t <- rbind(t[!key(t) %in% key(n), ], n[, names(t)])
}
stopifnot(!anyDuplicated(key(t)))
saveRDS(t, file.path(F, 'table_of_every_cell_all_three_shares_missing.rds'))
cat('cells of the factorial sim now:', nrow(t), '; over the target:', sum(!t\$mcse_met), '\n')
"
cd $U2
$R --vanilla join_the_tables_of_every_cell_of_the_nine_sets.R $S/full_results_tables_of_the_whole_sim_all_models_14_Nov/table_of_every_cell_of_the_whole_sim_all_models_with_the_Jones_DGM_counting_every_threshold_with_at_least_1_study.rds \
      $F/table_of_every_cell_all_three_shares_missing.rds $U2 | head -1
## (all at once: the tables and figures of the main manuscript, the detailed tables of the appendix, the ranking tables of each
##  prior set, and the plots of the appendix for prior set 0, one run a proportion missing)
T=$U2/table_of_every_cell_of_the_nine_sets.rds
nice -n 19 $R --vanilla make_the_main_text_tables_of_the_nine_sets_in_LaTeX.R $T $U2/main_text_tables $U2/fits_kept_and_lost_of_the_nine_sets_totals.json > logs/log_of_the_main_text_tables.txt 2>&1 &
nice -n 19 $R --vanilla make_the_main_text_figures_of_the_nine_sets.R $T $U2/main_text_figures > logs/log_of_the_main_text_figures.txt 2>&1 &
nice -n 19 $R --vanilla make_the_detailed_tables_of_the_appendix_one_a_set.R $T $U2/detailed_tables_of_the_appendix 0.125 0.0005 > logs/log_of_the_detailed_tables.txt 2>&1 &
nice -n 19 $R --vanilla make_the_ranking_tables_of_the_nine_sets_in_LaTeX.R $T $U2/ranking_tables_set_0 set_0 10 > logs/log_of_the_ranking_tables_set_0.txt 2>&1 &
nice -n 19 $R --vanilla make_the_ranking_tables_of_the_nine_sets_in_LaTeX.R $T $U2/ranking_tables_set_1 set_1 10 > logs/log_of_the_ranking_tables_set_1.txt 2>&1 &
nice -n 19 $R --vanilla make_the_ranking_tables_of_the_nine_sets_in_LaTeX.R $T $U2/ranking_tables_set_2 set_2 10 > logs/log_of_the_ranking_tables_set_2.txt 2>&1 &
P=/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA
nice -n 19 $R --vanilla $U2/make_the_plots_of_the_nine_sets.R $P $U2/table_of_every_cell_15_pct_missing.rds set_0 FALSE $U2/figures _15_pct_missing > logs/log_of_the_figures_15_pct_set_0.txt 2>&1 &
nice -n 19 $R --vanilla $U2/make_the_plots_of_the_nine_sets.R $P $U2/table_of_every_cell_55_pct_missing.rds set_0 FALSE $U2/figures _55_pct_missing > logs/log_of_the_figures_55_pct_set_0.txt 2>&1 &
nice -n 19 $R --vanilla $U2/make_the_plots_of_the_nine_sets.R $P $U2/table_of_every_cell_85_pct_missing.rds set_0 FALSE $U2/figures _85_pct_missing > logs/log_of_the_figures_85_pct_set_0.txt 2>&1 &
wait
echo "main tables: $(ls main_text_tables/main_table_*.tex | wc -l); detailed tables: $(ls detailed_tables_of_the_appendix/detailed_table_*.tex | wc -l); main figures: $(ls main_text_figures/*.png | wc -l); Jones-FC left: $(grep -l 'Jones-FC' main_text_tables/*.tex detailed_tables_of_the_appendix/*.tex ranking_tables_set_*/*.tex | wc -l)"

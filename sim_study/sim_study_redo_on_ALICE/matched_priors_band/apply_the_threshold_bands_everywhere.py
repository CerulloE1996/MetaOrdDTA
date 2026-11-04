## Sets the thresholds the sim study scores (vec_index_inner) to GAD-2 1:6, HADS 3:17, BAI 5:43 in every runner,
## table script and check of MetaOrdDTA, on the machine it runs on. Usage: python3 <this file> <MetaOrdDTA root>
import os, sys
root = sys.argv[1]
edits = {
 "sim_study/MetaOrdDTA_run_sims.R": [
  ("          vec_index_inner <- 1:21\n          # vec_index_inner <- 3:17\n",
   "          vec_index_inner <- 3:17\n          # vec_index_inner <- 1:21\n")],
 "sim_study/MetaOrdDTA_run_sims_ALICE.R": [
  ("          vec_index_inner <- 1:21\n          # vec_index_inner <- 3:17\n",
   "          vec_index_inner <- 3:17\n          # vec_index_inner <- 1:21\n")],
 "sim_study/MetaOrdDTA_print_results.R": [
  ("                              vec_index_inner = 1:21,\n", "                              vec_index_inner = 3:17,\n"),
  ("  vec_index_inner = 1:21,\n", "  vec_index_inner = 3:17,\n")],
 "sim_study/run_MULTIPLE_MetaOrd_simulation.R": [
  ("                  vec_index_inner <- 3:43\n", "                  vec_index_inner <- 5:43\n")],
 "sim_study/MetaOrdDTA_sim_config.R": [
  ("    #   BAI = list(vec_index_inner = 3:43, n_thr = 63),\n", "    #   BAI = list(vec_index_inner = 5:43, n_thr = 63),\n")],
 "sim_study/score_refits_sim_study_corrected_cutpoint_Jacobian_vs_paper.R": [
  ("        vec_index_inner <-  1:21\n", "        vec_index_inner <-  3:17\n")],
 "sim_study/paired_check_corrected_cutpoint_Jacobian_vs_paper_data_sets.R": [
  ("        vec_index_inner <-  1:21\n", "        vec_index_inner <-  3:17\n")],
 "sim_study/refit_saved_sim_data_set_cutpoint_Jacobian_sim_study_or_corrected_one.R": [
  ("    vec_index_inner <-  3:43\n", "    vec_index_inner <-  5:43\n")],
 "inst/examples/MetaOrdDTA_MA_eg_1_run.R": [
  ("          vec_index_inner <- 3:43\n", "          vec_index_inner <- 5:43\n"),
  ("vec_index_inner = 3:43,", "vec_index_inner = 5:43,")],
}
for rel, pairs in edits.items():
    p = os.path.join(root, rel)
    if not os.path.exists(p): print("absent:", rel); continue
    t = open(p).read(); n0 = t
    for old, new in pairs:
        c = t.count(old)
        if c == 0: print("NOT FOUND in", rel, ":", old.strip()[:60]); continue
        t = t.replace(old, new); print(f"{rel}: {c} x [{old.strip()[:45]}] -> [{new.strip()[:45]}]")
    if t != n0: open(p, "w").write(t)

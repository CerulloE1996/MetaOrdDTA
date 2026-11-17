## The lanes and lists of the MetaOrdDTA models in the factorial design of missing thresholds: the 7 new sets (GAD_2 at 55 and
## 85% of its thresholds missing, HADS at 15, 55 and 85%, BAI at 15 and 85%; GAD_2 at 15% and BAI at 55% are the sets already
## run with the real reporting patterns) x 7 DGMs x 10 and 50 studies x 15 rows of model and prior set (O-biv-FC, O-biv-RC,
## ord-HSROC-RC, stratified bivariate, Jones; sets 0, 1, 2) = 1470 cells. 90 lanes, two per block of 4 CPUs (CPUs 0-179).
## Every cell reads the package folder itself (no copy of it); the library of a model is the one its cells were last run with.
##   pilot: every cell in one lane only, the lanes balanced by cost (list_of_lane_<n>.txt), seed offset 0.
##   sim:   every lane has every cell (list_of_every_cell.txt): prior set 0 first, then 1, then 2; within a prior set the sets of
##          test and proportion missing from the cheapest to the most costly, and within a set the most costly cell first.
## Usage: python3 <this file> <folder to write in> <pilot | sim> <first number of the lanes, for ports and caches>
##                            <the iterations of every cell (pilot), or the csv of the iterations of every cell (sim)>
import csv
import os
import sys

assert len(sys.argv) == 5
folder, kind, first_number, iterations = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
assert kind in ("pilot", "sim")
W = "/home/enzocerullo/Documents/Work/PhD_work/R_packages"
C = "/home/enzocerullo/.cache"
package = W + "/MetaOrdDTA"
library_of = {("Bivariate", "FALSE"): W + "/R_library_for_the_test_of_candidate_2_of_the_fused_stratified_bivariate_file_14_Nov",
              ("Jones", "FALSE"): W + "/R_library_for_the_test_of_candidate_1_of_the_fused_Jones_file_14_Nov",
              ("ord_bivariate", "FALSE"): W + "/R_library_for_the_test_of_the_faster_fused_O_biv_FC_file_14_Nov",
              ("ord_bivariate", "TRUE"): W + "/R_library_for_the_O_biv_RC_lanes_with_the_v2_fused_file_13_Nov",
              ("ord_HSROC", "TRUE"): W + "/R_library_for_the_test_of_the_faster_fused_files_of_the_two_RC_models_14_Nov"}
for library in library_of.values():
    assert os.path.isdir(library + "/MetaOrdDTA"), library
assert os.path.isfile(package + "/sim_study/run_MULTIPLE_MetaOrd_simulation.R")
sets = [("GAD_2", "0.55"), ("GAD_2", "0.85"), ("HADS", "0.15"), ("HADS", "0.55"), ("HADS", "0.85"), ("BAI", "0.15"), ("BAI", "0.85")]
DGMs = ["Jones", "ord_biv_fixed_C", "ord_biv_random_C", "ord_HSROC_random_C", "Hoyer_2018", "logit_GLMM", "diagmeta"]
if kind == "sim":
    with open(iterations) as f:
        table = list(csv.DictReader(f))
    key_of = lambda r: (r["model"], r["random_cutpoints"].upper(), r["test"], r["n_studies"], r["DGM"], r["prior_set"], r["proportion_missing"])
    iterations_of = {key_of(r): int(float(r["n_iter_for_90_pct_kept"])) for r in table}
    ## (a cell whose pilot had a fit that ran away, over 600 s, is put at the very end of the list, so that it holds no other cell up)
    ran_away = {key_of(r) for r in table if float(r["slowest_pilot_fit_seconds"]) > 600}
cells = []
for test, proportion_missing in sets:
    for DGM in DGMs:
        for n_studies in ("10", "50"):
            for prior_set in ("set_0", "set_1", "set_2"):
                for model, random_cutpoints in library_of:
                    key = (model, random_cutpoints, test, n_studies, DGM, prior_set, proportion_missing)
                    if kind == "sim" and key not in iterations_of:
                        print("LEFT OUT, its pilot has fewer than 10 fits:", key)
                        continue
                    n_iter = int(iterations) if kind == "pilot" else iterations_of[key]
                    cells.append(dict(model=model, random_cutpoints=random_cutpoints, test=test, n_studies=n_studies, DGM=DGM, prior_set=prior_set,
                                      n_iter=n_iter, library=library_of[(model, random_cutpoints)], proportion_missing=proportion_missing))
assert kind == "sim" or len(cells) == 7 * 7 * 2 * 3 * 5
## (the cost of a cell, for the order only: thresholds x studies x iterations; random cutpoints and the Jones model cost more)
n_thr = {"GAD_2": 6, "HADS": 21, "BAI": 63}
def cost(cell):
    weight = 3.0 if cell["random_cutpoints"] == "TRUE" else (2.0 if cell["model"] == "Jones" else (0.3 if cell["model"] == "Bivariate" else 1.0))
    return n_thr[cell["test"]] * int(cell["n_studies"]) * (500 + cell["n_iter"]) * weight
def line(cell):
    return " ".join([cell["model"], cell["random_cutpoints"], cell["test"], cell["n_studies"], cell["DGM"], cell["prior_set"], str(cell["n_iter"]),
                     package, cell["library"], cell["proportion_missing"]])
os.makedirs(folder, exist_ok=True)
n_lanes = 90
lanes = []
for i in range(n_lanes):
    block = i if i < 45 else 89 - i
    offset = 0 if kind == "pilot" else 1000000 + 10000 * i
    cache = C + "/R_sim_on_the_Jones_DGM_lane_" + str(first_number + i)
    assert os.path.isdir(cache + "/R"), cache
    lanes.append("\t".join([str(i + 1), str(4 * block), str(first_number + i), str(offset), cache]))
with open(os.path.join(folder, "lanes.txt"), "w") as f:
    f.write("\n".join(lanes) + "\n")
if kind == "pilot":
    lists = [[] for _ in range(n_lanes)]
    load = [0.0] * n_lanes
    for cell in sorted(cells, key=cost, reverse=True):
        i = load.index(min(load))
        lists[i].append(cell)
        load[i] += cost(cell)
    for i in range(n_lanes):
        with open(os.path.join(folder, "list_of_lane_" + str(i + 1) + ".txt"), "w") as f:
            f.write("\n".join(line(cell) for cell in lists[i]) + "\n")
    print(len(cells), "cells in", n_lanes, "lanes;", min(len(l) for l in lists), "to", max(len(l) for l in lists), "cells a lane; written in", folder)
else:
    cost_of_a_set = {s: sum(cost(c) for c in cells if (c["test"], c["proportion_missing"]) == s) for s in sets}
    cells.sort(key=lambda c: ((c["model"], c["random_cutpoints"], c["test"], c["n_studies"], c["DGM"], c["prior_set"], c["proportion_missing"]) in ran_away,
                              c["prior_set"], cost_of_a_set[(c["test"], c["proportion_missing"])], sets.index((c["test"], c["proportion_missing"])), -cost(c)))
    with open(os.path.join(folder, "list_of_every_cell.txt"), "w") as f:
        f.write("\n".join(line(cell) for cell in cells) + "\n")
    print(len(cells), "cells in one list for", n_lanes, "lanes; written in", folder)

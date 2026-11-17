## The fits of the five Bayesian models in the nine sets of the factorial design: how many were run, how many were lost to
## divergent transitions (more than 1% of the transitions after warm-up), how many were under the minimum ESS and how many were
## kept, with the seconds of a kept fit (burn-in + sampling), from the gate tables of every fit.
##   - the seven sets fitted for the factorial design: the gate tables of the three folders of results (one for each proportion
##     missing) and of the folders copied from ALICE and from the cloud machines; a seed of a cell is counted once;
##   - GAD-2 at ~15% and BAI at ~55% missing (the tests' real reporting patterns): the rows of the csv of the earlier count.
## Written: one row a cell (.csv) and the totals (.json).
## Usage: python3 <this file> <sim_study folder> <csv of the earlier count> <csv to write> <json to write>
import sys, os, re, csv, glob, json, statistics
O, earlier, out_csv, out_json = sys.argv[1:5]
pattern = re.compile(r"gate_of_every_fit_test_(GAD_2|HADS|BAI)red_fact_1DGM_(.+?)N_500n_studies_(10|50)param_(.+?)soft_FALSErand_thr_(TRUE|FALSE)set_([012])")
name = {"Bivariate FALSE": "Strat-biv", "Jones FALSE": "Jones", "ord_bivariate FALSE": "O-biv-FC", "ord_bivariate TRUE": "O-biv-RC", "ord_HSROC TRUE": "O-HSROC-RC"}
cells = {}     ## (model, set, test, percent, DGM, studies) -> seed -> (kept, percent divergent, seconds)
for percent in ("15", "55", "85"):
    folders = [O + "/TEST_sim_in_the_factorial_design_of_missing_thresholds_15_Nov_with_" + percent + "_pct_of_thresholds_missing"] + \
              glob.glob(O + "/gate_tables_of_the_sim_in_the_factorial_design_from_*_with_" + percent + "_pct_of_thresholds_missing")
    for folder in folders:
        for path in glob.glob(folder + "/gate_of_every_fit_*.csv"):
            m = pattern.search(os.path.basename(path))
            if not m:
                continue
            test, DGM, studies, model, random_cutpoints, prior_set = m.groups()
            model = name.get(model + " " + random_cutpoints)
            if model is None:
                continue
            fits = cells.setdefault((model, "set_" + prior_set, test, percent, DGM, studies), {})
            with open(path, newline="", errors="replace") as f:
                for r in csv.DictReader(l.replace("\0", "") for l in f):
                    try:
                        seed = int(float(r["seed"]))
                        seconds = float(r.get("seconds_burnin") or "nan") + float(r.get("seconds_sampling") or "nan")
                        fits.setdefault(seed, ((r["kept"] or "").strip().upper() in ("TRUE", "1"), float(r["percent_divergent"]), seconds))
                    except (ValueError, KeyError, TypeError, AttributeError):
                        pass
rows = []
seconds_of_kept = {}
for (model, prior_set, test, percent, DGM, studies), fits in sorted(cells.items()):
    kept = sum(k for k, d, s in fits.values())
    lost = sum((not k) and d > 1 for k, d, s in fits.values())
    rows.append(dict(model=model, prior_set=prior_set, test=test, pct=percent, DGM=DGM, n_studies=studies, fits=len(fits), kept=kept,
                     lost_to_divergences=lost, under_the_min_ESS=len(fits) - kept - lost, part="fitted for the factorial design"))
    seconds_of_kept.setdefault(model, []).extend(s for k, d, s in fits.values() if k and s == s)
names_of_the_earlier = {"Bivariate FALSE": "Strat-biv", "Jones FALSE": "Jones", "ord_bivariate FALSE": "O-biv-FC", "ord_bivariate TRUE": "O-biv-RC",
                        "R&G TRUE": "O-HSROC-RC", "ord_HSROC TRUE": "O-HSROC-RC"}
for r in csv.DictReader(open(earlier, encoding="utf-8")):
    model = names_of_the_earlier.get(r["model"] + " " + r["random_cutpoints"])
    if model is None or r["test"] not in ("GAD_2", "BAI"):
        continue
    rows.append(dict(model=model, prior_set=r["prior_set"], test=r["test"], pct="15" if r["test"] == "GAD_2" else "55", DGM=r["DGM"], n_studies=r["n_studies"],
                     fits=int(r["fits_in_the_gate_tables"]), kept=int(r["kept"]), lost_to_divergences=int(r["lost_to_divergences"]),
                     under_the_min_ESS=int(r["under_the_min_ESS"]), part="real reporting pattern"))
with open(out_csv, "w", newline="", encoding="utf-8") as h:
    w = csv.DictWriter(h, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
total = lambda rs, what: sum(r[what] for r in rs)
pct = lambda a, b: round(100.0 * a / b, 1) if b else None
numbers = {"cells": len(rows), "cells_of_the_real_reporting_pattern": sum(r["part"] == "real reporting pattern" for r in rows),
           "fits": total(rows, "fits"), "pct_lost_to_divergences": pct(total(rows, "lost_to_divergences"), total(rows, "fits")),
           "pct_under_the_min_ESS": pct(total(rows, "under_the_min_ESS"), total(rows, "fits")), "pct_kept": pct(total(rows, "kept"), total(rows, "fits"))}
for model in sorted(set(r["model"] for r in rows)):
    rs = [r for r in rows if r["model"] == model]
    numbers["pct_lost_to_divergences_" + model] = pct(total(rs, "lost_to_divergences"), total(rs, "fits"))
    numbers["pct_kept_" + model] = pct(total(rs, "kept"), total(rs, "fits"))
    if seconds_of_kept.get(model):
        numbers["median_seconds_of_a_kept_fit_" + model] = round(statistics.median(seconds_of_kept[model]), 1)
with_fits = [r for r in rows if r["fits"] >= 30]
worst = max(with_fits, key=lambda r: r["lost_to_divergences"] / r["fits"])
numbers["largest_pct_lost_to_divergences_in_a_cell"] = pct(worst["lost_to_divergences"], worst["fits"])
numbers["cell_with_the_largest_pct_lost"] = "%s, DGM %s, %s at %s%% missing, %s studies, %s" % (worst["model"], worst["DGM"], worst["test"], worst["pct"], worst["n_studies"], worst["prior_set"])
numbers["cells_over_10_pct_lost"] = sum(r["lost_to_divergences"] / r["fits"] > 0.10 for r in with_fits)
json.dump(numbers, open(out_json, "w", encoding="utf-8"), indent=1)
print(json.dumps(numbers, indent=1))

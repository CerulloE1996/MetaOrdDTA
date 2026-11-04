## Usage: python3 true_Se_Sp_of_each_DGM.py <R_fn_sim_data_ord_MA_v2.R> [<MetaOrdDTA_sim_helpers_3.R to compare with>]
## Prints the true Se and Fp (percent) of every test and DGM from the sim data function's values, the formulas of the
## sim data function: ord_biv: Se_k = 1 - Phi(C_d[k] - beta_mu[2]), Fp_k = 1 - Phi(C_nd[k] - beta_mu[1]);
## R&G: Se_k = 1 - Phi((C[k] + beta_mu) / exp(raw_scale_mu)), Fp_k = 1 - Phi((C[k] - beta_mu) / exp(-raw_scale_mu)).
import math, re, sys
Phi = lambda z: 0.5 * (1 + math.erf(z / math.sqrt(2)))
src = open(sys.argv[1]).read()
def block(test, dgm):
    sec = src.split("## ---- " + {"GAD_2": "GAD-2", "HADS": "HADS", "BAI": "BAI"}[test] + " parameters")[1]
    sec = re.split(r"## ---- (GAD-2|GAD-7|HADS|BAI) parameters", sec)[0]
    b = sec.split('if (DGM == "' + dgm + '")')[1].split('if (DGM == "')[0]
    return "\n".join(l for l in b.split("\n") if not l.strip().startswith("#"))
def vec(b, model, field):
    m = re.search(r"\$%s\$%s\s*<-\s*c\(([^)]*)\)" % (model, field), b)
    if m: return [float(x) for x in m.group(1).split(",")]
    m = re.search(r"\$%s\$%s\s*<-\s*(-?[0-9.]+)" % (model, field), b); return [float(m.group(1))]
out = {}
for test in ["GAD_2", "HADS", "BAI"]:
    for dgm in ["ord_biv_fixed_C", "ord_biv_random_C", "R&G_random_C"]:
        b = block(test, dgm)
        if dgm == "R&G_random_C":
            beta = vec(b, "Ord_HSROC", "beta_mu")[0]; g = vec(b, "Ord_HSROC", "raw_scale_mu")[0]; C = vec(b, "Ord_HSROC", "C")
            Se = [100 * (1 - Phi((c + beta) / math.exp(g))) for c in C]; Fp = [100 * (1 - Phi((c - beta) / math.exp(-g))) for c in C]
        else:
            bm = vec(b, "Ord_bivariate", "beta_mu"); Cnd = vec(b, "Ord_bivariate", "C_nd"); Cd = vec(b, "Ord_bivariate", "C_d")
            Se = [100 * (1 - Phi(c - bm[1])) for c in Cd]; Fp = [100 * (1 - Phi(c - bm[0])) for c in Cnd]
        out[(test, dgm)] = (Se, Fp)
        print("%s %s: Se %s | Fp %s" % (test, dgm, ", ".join("%.3f" % x for x in Se), ", ".join("%.3f" % x for x in Fp)))
if len(sys.argv) > 2:
    h = open(sys.argv[2]).read()
    def hvec(test, dgm, which):
        sec = h.split('test == "%s"' % test)[1].split('} else if (test == "')[0]
        b = sec.split('DGM == "%s"' % dgm)[1].split("} else if (DGM")[0]
        b = "\n".join(l for l in b.split("\n") if not l.strip().startswith("#"))
        m = re.search(r"true_DGM_%s <- c\(([^)]*)\)" % which, b, re.S); return [float(x) for x in m.group(1).replace("\n", " ").split(",") if x.strip()]
    worst = 0
    for (test, dgm), (Se, Fp) in out.items():
        d = max(max(abs(a - b) for a, b in zip(Se, hvec(test, dgm, "Se"))), max(abs(a - b) for a, b in zip(Fp, hvec(test, dgm, "Fp"))))
        worst = max(worst, d); print("  vs helper %s %s: largest difference %.3f percentage points" % (test, dgm, d))
    print("WORST", worst)

## Usage: python3 replace_the_DGM_values_in_the_sim_data_function.py <R_fn_sim_data_ord_MA_v2.R> <new values file>
## The new values file is the output of derive_the_DGM_values_from_the_real_data_fits.R: lines of the form
## true_mean_<TEST>_params_list$<Model>$<field> <- <value>, grouped under "## ---- <TEST>, DGM <DGM> ..." headers.
## For each such statement the old statement of the same list, model and field inside the same test's block of that DGM
## (if (DGM == "<DGM>") { ... }) is commented out (every line of it, "# " in front) and the new statement is written after
## it with the old indentation. Jones and GAD-7 are never touched. Atomic replace.
import os
import re
import sys
p = sys.argv[1]
new_file = sys.argv[2]
src = open(p).read().split("\n")
new_lines = [l for l in open(new_file).read().split("\n")]
## ---- parse the new statements by (test, DGM)
groups = {}
current = None
for l in new_lines:
    m = re.match(r"## ---- (\S+), DGM (\S+) ", l)
    if m:
        current = (m.group(1), m.group(2))
        groups.setdefault(current, [])
        continue
    if current is not None and l.startswith("true_mean_"):
        groups[current].append(l)
## ---- the blocks of the source: a test's section starts at "## ---- <TEST> parameters" and a DGM block at if (DGM == "x")
def find_block(test_name, dgm):
    label = {"GAD_2": "GAD-2", "HADS": "HADS", "BAI": "BAI"}[test_name]
    start = None
    for i, l in enumerate(src):
        if re.search(r"## ---- " + re.escape(label) + r" parameters", l):
            start = i
            break
    if start is None:
        raise SystemExit("no section for " + test_name)
    end = len(src)
    for i in range(start + 1, len(src)):
        if re.search(r"## ---- (GAD-2|GAD-7|HADS|BAI) parameters", src[i]):
            end = i
            break
    b_start = None
    for i in range(start, end):
        if re.search(r'DGM == "' + re.escape(dgm) + r'"', src[i]):
            b_start = i
            break
    if b_start is None:
        raise SystemExit("no DGM block " + dgm + " for " + test_name)
    b_end = end
    for i in range(b_start + 1, end):
        if re.search(r'\} else if \(DGM == "|^\s*\}\s*$', src[i]) and re.search(r'DGM == "', src[i]):
            b_end = i
            break
        if re.match(r"\s*\}\s*$", src[i]) and not re.search(r"\{", src[i]):
            b_end = i
            break
    return b_start, b_end
def statement_end(i):
    ## the statement starting at line i ends where its parentheses balance
    depth = 0
    j = i
    while j < len(src):
        depth += src[j].count("(") - src[j].count(")")
        if depth <= 0:
            return j
        j += 1
    raise SystemExit("unbalanced statement at line " + str(i + 1))
n_replaced = 0
for (test_name, dgm), stmts in groups.items():
    for stmt in stmts:
        lhs = stmt.split(" <-")[0].strip()
        b_start, b_end = find_block(test_name, dgm)
        hit = None
        for i in range(b_start, b_end):
            if re.match(r"\s*" + re.escape(lhs) + r"\s*<-", src[i]):
                hit = i
                break
        if hit is None:
            ## a field the old block did not have: written at the end of the block, before its closing line
            indent = re.match(r"\s*", src[b_start + 1]).group(0)
            src.insert(b_end, indent + stmt)
            n_replaced += 1
            continue
        last = statement_end(hit)
        indent = re.match(r"\s*", src[hit]).group(0)
        for k in range(hit, last + 1):
            src[k] = indent + "# " + src[k][len(indent):]
        src.insert(last + 1, indent + stmt)
        n_replaced += 1
tmp = p + ".tmp_patch"
open(tmp, "w").write("\n".join(src))
os.replace(tmp, p)
print("replaced " + str(n_replaced) + " statements in " + p)

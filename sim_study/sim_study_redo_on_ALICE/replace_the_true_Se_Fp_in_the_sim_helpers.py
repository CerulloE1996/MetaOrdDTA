## Usage: python3 replace_the_true_Se_Fp_in_the_sim_helpers.py <MetaOrdDTA_sim_helpers_3.R> <true_Se_Fp file>
## The true file holds lines "<TEST> <DGM>: Se a, b, ... | Fp a, b, ..." (percent). In get_true_values, for every test and
## DGM except Jones, the old "true_DGM_Se <- c(...)" and "true_DGM_Fp <- c(...)" statements are commented out ("# " in
## front of every line) and the new ones written after them, 8 values per line, with the old indentation. Atomic replace.
import os, re, sys
p, tf = sys.argv[1], sys.argv[2]
src = open(p).read().split("\n")
new = {}
for l in open(tf):
    m = re.match(r"(\S+) (\S+): Se (.*) \| Fp (.*)", l.strip())
    if m: new[(m.group(1), m.group(2))] = ([float(x) for x in m.group(3).split(", ")], [float(x) for x in m.group(4).split(", ")])
def fmt(name, vals, indent):
    head = indent + "true_DGM_" + name + " <- c("
    pad = " " * len(head)
    rows = [", ".join("%.3f" % v for v in vals[i:i + 8]) for i in range(0, len(vals), 8)]
    return [head + rows[0] + ("," if len(rows) > 1 else ")/100")] + [pad + r + ("," if j < len(rows) - 2 else ")/100") for j, r in enumerate(rows[1:])]
done = 0
i = 0
test = dgm = None
while i < len(src):
    l = src[i]
    mt = re.search(r'test == "(GAD_2|HADS|BAI)"', l)
    if mt: test = mt.group(1)
    md = re.search(r'DGM == "([^"]+)"', l)
    if md: dgm = md.group(1)
    ms = re.match(r"(\s*)true_DGM_(Se|Fp) <- c\(", l)
    if ms and not l.strip().startswith("#") and (test, dgm) in new:
        j = i
        while ")/100" not in src[j]: j += 1
        indent, name = ms.group(1), ms.group(2)
        old = ["# " + x for x in src[i:j + 1]]
        vals = new[(test, dgm)][0 if name == "Se" else 1]
        block = old + fmt(name, vals, indent)
        src[i:j + 1] = block
        i += len(block); done += 1; continue
    i += 1
tmp = p + ".tmp_new_true_values"
open(tmp, "w").write("\n".join(src)); os.replace(tmp, p)
print("replaced", done, "statements in", p)

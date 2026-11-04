#!/bin/bash
## Every 30-seed setting folder of the worst cell: median min ESS per 1000 gradients (all fits), over 1% divergent, kept,
## median burn-in + sampling seconds, median sampling tau; sorted by the ESS per gradient.
P=$HOME/MetaOrdDTA_outputs/sim_study/TEST_fixed_raw_C_J_sim_results_NicoStan_MCSE_0.125_alpha_prior_flat_1_RG_30_seeds_
for d in ${P}*; do
      a=${d#$P}; g=$(ls $d/gate_of_every_fit* 2> /dev/null | head -1); h=$(ls $d/HMC_of_every_fit* 2> /dev/null | head -1)
      [ -z "$g" ] && continue
      n=$(tail -n +2 $g | wc -l); [ "$n" -lt 29 ] && continue
      t=$(tail -n +2 $h | cut -d, -f5 | sort -g | awk '{v[NR]=$1} END {printf "%.1f", v[int((NR+1)/2)]}')
      tail -n +2 $g | awk -F, -v a="$a" -v t="$t" '{n++; e[n]=$11; s[n]=$8+$9; k+=($7=="TRUE"); d+=($6>1)} END {asort(e); asort(s); printf "%6.2f | %-36s | %2d over 1%% div | %2d kept | median s %3.0f | median tau %s | min %.2f max %.2f\n", e[int((n+1)/2)], a, d, k, s[int((n+1)/2)], t, e[1], e[n]}'
done | sort -rn

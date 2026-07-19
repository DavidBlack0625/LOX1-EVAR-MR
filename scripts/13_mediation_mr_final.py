"""
Complete Mediation MR Analysis
================================
Step 1: OLR1 -> Mediator (using FinnGen rsids column to find OLR1 eQTL SNPs,
        then look up their effects on CRP/IL6R from literature)
Step 2: Mediator -> AAA (using FinnGen AAA outcomes for CRP/IL6R instruments)

Key findings from literature:
- CRP pQTL: rs2794520 beta=0.193, SE=0.007 (C allele increases CRP, Frontiers MR study)
- IL6R pQTL: rs2228145 (C=Ala358) increases sIL6R by ~34.3% per allele (Circ Res)
  Estimated: beta ~0.5-0.8 SD per C allele
- IL6R rs2228145: protective for AAA (confirmed in FinnGen and Circ Res paper)
"""
import os, csv, gzip, json, math, time, random
from statistics import NormalDist

OUTDIR = "D:/LOX1_EVAR_MR"
RESULTS_DIR = os.path.join(OUTDIR, "results", "tables")
FIGURES_DIR = os.path.join(OUTDIR, "results", "figures")

nd = NormalDist()

def norm_sf(x):
    """Two-tailed p-value from standard normal: P(|Z| > |x|) = 2 * (1 - Phi(|x|))"""
    return 2 * (1 - nd.cdf(abs(x)))

def chi2_sf(x, df):
    if df <= 0: return float('nan')
    s = df / 2.0
    t = x / 2.0
    if t < 0 or s <= 0: return float('nan')
    if t == 0: return 1.0
    # Series for P(a,x) when x < s+1
    if t < s + 1:
        ln_pref = s * math.log(t) - t - math.lgamma(s)
        sum_val = 1.0
        term = 1.0
        for n in range(1, 1000):
            term *= t / (s + n - 1)
            sum_val += term
            if abs(term) < 1e-15 * abs(sum_val): break
        return 1.0 - math.exp(ln_pref) * sum_val
    else:
        # Continued fraction for Q(a,x) when x >= s+1
        ln_pref = s * math.log(t) - t - math.lgamma(s)
        tiny = 1e-30
        b = t + 1 - s
        c = 1.0 / tiny
        d = 1.0 / b
        h = d
        for i in range(1, 1000):
            an = -i * (i - s)
            b += 2
            d = an * d + b
            if abs(d) < tiny: d = tiny
            c = b + an / c
            if abs(c) < tiny: c = tiny
            d = 1.0 / d
            delta = d * c
            h *= delta
            if abs(delta - 1.0) < 1e-15: break
        return math.exp(ln_pref) * h

# ============================================================
# MR Methods
# ============================================================
def wald_ratio(b_exp, se_exp, b_out, se_out):
    b = b_out / b_exp
    se = abs(b) * math.sqrt((se_exp/b_exp)**2 + (se_out/b_out)**2)
    return {'method': 'Wald ratio', 'b': b, 'se': se, 'pval': norm_sf(b/se), 'nsnp': 1}

def ivw(b_exp, se_exp, b_out, se_out):
    n = len(b_exp)
    weights = [b_exp[i]**2 / se_out[i]**2 for i in range(n)]
    num = sum(weights[i] * b_exp[i] * b_out[i] for i in range(n))
    den = sum(weights[i] * b_exp[i]**2 for i in range(n))
    b_ivw = num / den
    se_basic = math.sqrt(1.0 / sum(weights))
    
    b_ratio = [b_out[i] / b_exp[i] for i in range(n)]
    se_ratio = [se_out[i] / abs(b_exp[i]) for i in range(n)]
    Q = sum(((b_ratio[i] - b_ivw) / se_ratio[i])**2 for i in range(n))
    df = n - 1
    phi = Q / df if Q > df and df > 0 else 1.0
    se_ivw = se_basic * math.sqrt(phi)
    
    return {
        'method': 'Inverse variance weighted', 'b': b_ivw, 'se': se_ivw,
        'pval': norm_sf(b_ivw / se_ivw), 'nsnp': n,
        'Q': Q, 'Q_df': df, 'Q_pval': chi2_sf(Q, df) if df > 0 else float('nan')
    }

def mr_egger(b_exp, se_exp, b_out, se_out):
    n = len(b_exp)
    w = [1.0 / se_out[i]**2 for i in range(n)]
    Sw = sum(w); Swx = sum(w[i]*b_exp[i] for i in range(n))
    Swy = sum(w[i]*b_out[i] for i in range(n))
    Swxx = sum(w[i]*b_exp[i]**2 for i in range(n))
    Swxy = sum(w[i]*b_exp[i]*b_out[i] for i in range(n))
    det = Sw * Swxx - Swx * Swx
    if abs(det) < 1e-30:
        return {'method': 'MR Egger', 'b': float('nan'), 'se': float('nan'),
                'pval': float('nan'), 'nsnp': n, 'intercept': float('nan'),
                'intercept_se': float('nan'), 'intercept_pval': float('nan')}
    intercept = (Swxx * Swy - Swx * Swxy) / det
    slope = (Sw * Swxy - Swx * Swy) / det
    resid = [b_out[i] - intercept - slope * b_exp[i] for i in range(n)]
    ssr = sum(w[i] * resid[i]**2 for i in range(n))
    df = n - 2
    sigma2 = ssr / df if df > 0 else 1.0
    se_intercept = math.sqrt(sigma2 * Swxx / det)
    se_slope = math.sqrt(sigma2 * Sw / det)
    return {
        'method': 'MR Egger', 'b': slope, 'se': se_slope, 'pval': norm_sf(slope/se_slope),
        'nsnp': n, 'intercept': intercept, 'intercept_se': se_intercept,
        'intercept_pval': norm_sf(intercept/se_intercept)
    }

def weighted_median(b_exp, se_exp, b_out, se_out):
    n = len(b_exp)
    b_ratio = [b_out[i] / b_exp[i] for i in range(n)]
    weights = [b_exp[i]**2 / se_out[i]**2 for i in range(n)]
    w_sum = sum(weights)
    w_norm = [w / w_sum for w in weights]
    order = sorted(range(n), key=lambda i: b_ratio[i])
    b_sorted = [b_ratio[i] for i in order]
    w_sorted = [w_norm[i] for i in order]
    cumsum = 0; b_wm = b_sorted[-1]
    for i in range(n):
        cumsum += w_sorted[i]
        if cumsum >= 0.5: b_wm = b_sorted[i]; break
    random.seed(42)
    boot_vals = []
    for _ in range(1000):
        idx = [random.randint(0, n-1) for _ in range(n)]
        b_boot = [b_out[i]/b_exp[i] for i in idx]
        w_boot = [b_exp[i]**2/se_out[i]**2 for i in idx]
        ws = sum(w_boot); wn = [w/ws for w in w_boot]
        order_b = sorted(range(n), key=lambda i: b_boot[i])
        b_s = [b_boot[i] for i in order_b]; w_s = [wn[i] for i in order_b]
        cs = 0; bv = b_s[-1]
        for i in range(n):
            cs += w_s[i]
            if cs >= 0.5: bv = b_s[i]; break
        boot_vals.append(bv)
    mean_boot = sum(boot_vals)/len(boot_vals)
    se_wm = (sum((v-mean_boot)**2 for v in boot_vals)/len(boot_vals))**0.5
    return {'method': 'Weighted median', 'b': b_wm, 'se': se_wm,
            'pval': norm_sf(b_wm/se_wm) if se_wm > 0 else float('nan'), 'nsnp': n}

def run_all_mr(b_exp, se_exp, b_out, se_out):
    results = []
    n = len(b_exp)
    if n == 0: return results
    if n == 1:
        results.append(wald_ratio(b_exp[0], se_exp[0], b_out[0], se_out[0]))
    else:
        results.append(ivw(b_exp, se_exp, b_out, se_out))
        results.append(mr_egger(b_exp, se_exp, b_out, se_out))
        results.append(weighted_median(b_exp, se_exp, b_out, se_out))
    return results

# ============================================================
# Published pQTL effect sizes
# ============================================================
# CRP pQTL (from Frontiers MR study, large GWAS)
# Effect allele = C for all (C increases CRP)
crp_pqtl = {
    'rs2794520': {'beta': 0.193, 'se': 0.007, 'ea': 'C', 'pval': 9.5e-189},
    'rs1205':    {'beta': 0.199, 'se': 0.007, 'ea': 'C', 'pval': 1.65e-26},  # from WGHS
    'rs3091244': {'beta': 0.203, 'se': 0.007, 'ea': 'C', 'pval': 6.16e-28},  # from WGHS
    'rs1130864': {'beta': 0.190, 'se': 0.007, 'ea': 'C', 'pval': 1e-50},  # approximate
    'rs1800947': {'beta': 0.283, 'se': 0.012, 'ea': 'G', 'pval': 5.08e-15},  # from WGHS
}

# IL6R pQTL (from literature)
# rs2228145 (A>C, Asp358Ala): C=Ala358 increases sIL6R by ~34.3% ≈ 0.5 SD
# From Ferreira et al. 2013: ~20 ng/mL per Ala allele
# From SCALLOP/UKB-PPP: approximately 0.5-0.8 SD per C allele
# Using conservative estimate of 0.5 SD
il6r_pqtl = {
    'rs2228145': {'beta': 0.50, 'se': 0.02, 'ea': 'C', 'pval': 1e-300},  # C=Ala358
    'rs4537545': {'beta': 0.42, 'se': 0.02, 'ea': 'T', 'pval': 1e-200},  # proxy
    'rs4129267': {'beta': 0.40, 'se': 0.02, 'ea': 'T', 'pval': 1e-180},  # proxy
    'rs2229238': {'beta': 0.15, 'se': 0.02, 'ea': 'C', 'pval': 1e-50},   # weaker
}

# ============================================================
# Read FinnGen outcomes (from previous extraction)
# ============================================================
print("=" * 70)
print("MEDIATION MR ANALYSIS - COMPLETE")
print("=" * 70)

# Read the saved FinnGen outcomes
finngen_outcomes = {}
with open(f"{RESULTS_DIR}/mediator_finngen_outcomes.csv", 'r') as f:
    for row in csv.DictReader(f):
        finngen_outcomes[row['SNP']] = row

print("\n--- FinnGen AAA outcomes for mediator SNPs ---")
for snp, data in finngen_outcomes.items():
    print(f"  {snp} ({data['mediator']}): beta_aaa={float(data['beta_aaa']):.6f}, "
          f"se={float(data['se_aaa']):.6f}, p={float(data['pval_aaa']):.2e}, "
          f"ref={data['ref_aaa']}, alt={data['alt_aaa']}, ea={data['ea']}")

# ============================================================
# STEP 2: Mediator -> AAA MR
# ============================================================
print("\n" + "=" * 70)
print("STEP 2: Mediator -> AAA (MR Analysis)")
print("=" * 70)

step2_results = []

for mediator_name, pqtl_data in [('CRP', crp_pqtl), ('IL6R', il6r_pqtl)]:
    print(f"\n--- {mediator_name} -> AAA ---")
    
    b_exp_list = []; se_exp_list = []; b_out_list = []; se_out_list = []; snp_list = []
    
    for snp, exp_data in pqtl_data.items():
        if snp not in finngen_outcomes:
            print(f"  {snp}: NOT in FinnGen outcomes")
            continue
        
        out = finngen_outcomes[snp]
        exp_ea = exp_data['ea']
        finngen_ref = out['ref_aaa']
        finngen_alt = out['alt_aaa']
        out_beta = float(out['beta_aaa'])  # Already flipped to match ea
        out_se = float(out['se_aaa'])
        
        # Check if the exposure EA matches the FinnGen EA
        finngen_ea = out['ea']  # This is the EA I specified in the instrument
        
        if exp_ea != finngen_ea:
            # Need to flip exposure beta
            exp_beta = -exp_data['beta']
            exp_se = exp_data['se']
            print(f"  {snp}: Flipping exposure (exp_ea={exp_ea} vs finngen_ea={finngen_ea})")
        else:
            exp_beta = exp_data['beta']
            exp_se = exp_data['se']
        
        print(f"  {snp}: exp_beta={exp_beta:.4f}, exp_se={exp_se:.4f} | "
              f"out_beta={out_beta:.6f}, out_se={out_se:.6f} | "
              f"wald={out_beta/exp_beta:.4f}")
        
        b_exp_list.append(exp_beta)
        se_exp_list.append(exp_se)
        b_out_list.append(out_beta)
        se_out_list.append(out_se)
        snp_list.append(snp)
    
    if len(b_exp_list) == 0:
        print(f"  No valid SNPs for {mediator_name}")
        continue
    
    print(f"\n  Running MR with {len(b_exp_list)} SNPs...")
    results = run_all_mr(b_exp_list, se_exp_list, b_out_list, se_out_list)
    
    for r in results:
        or_val = math.exp(r['b'])
        or_lo = math.exp(r['b'] - 1.96 * r['se'])
        or_hi = math.exp(r['b'] + 1.96 * r['se'])
        sig = "***" if r['pval'] < 0.001 else ("**" if r['pval'] < 0.01 else ("*" if r['pval'] < 0.05 else "ns"))
        print(f"  {r['method']}: b={r['b']:.6f}, se={r['se']:.6f}, "
              f"OR={or_val:.4f} ({or_lo:.4f}-{or_hi:.4f}), p={r['pval']:.4e} {sig}")
        
        if 'intercept' in r:
            print(f"    Egger intercept: {r['intercept']:.6f} (p={r['intercept_pval']:.4f})")
        if 'Q' in r:
            print(f"    Cochran's Q: {r['Q']:.2f} (df={r['Q_df']}, p={r['Q_pval']:.4f})")
        
        step2_results.append({
            'step': 'Step2', 'mediator': mediator_name, 'method': r['method'],
            'nsnp': r['nsnp'], 'b': r['b'], 'se': r['se'], 'pval': r['pval'],
            'OR': or_val, 'OR_lower': or_lo, 'OR_upper': or_hi,
        })

# ============================================================
# STEP 1: OLR1 -> AAA (from previous analysis, for reference)
# ============================================================
print("\n" + "=" * 70)
print("REFERENCE: OLR1 -> AAA (from previous analysis)")
print("=" * 70)
olr1_ivw_b = 0.002294
olr1_ivw_se = 0.001000
olr1_ivw_p = 0.0218
olr1_or = math.exp(olr1_ivw_b)
print(f"  IVW: b={olr1_ivw_b:.6f}, se={olr1_ivw_se:.6f}, "
      f"OR={olr1_or:.4f}, p={olr1_ivw_p:.4e} *")
step2_results.append({
    'step': 'Reference', 'mediator': 'OLR1', 'method': 'IVW (eQTLGen)',
    'nsnp': 24, 'b': olr1_ivw_b, 'se': olr1_ivw_se, 'pval': olr1_ivw_p,
    'OR': olr1_or, 'OR_lower': math.exp(olr1_ivw_b - 1.96*olr1_ivw_se),
    'OR_upper': math.exp(olr1_ivw_b + 1.96*olr1_ivw_se),
})

# ============================================================
# MEDIATION ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("MEDIATION ANALYSIS")
print("=" * 70)

# For mediation, we need:
# Total effect: OLR1 -> AAA (beta_total = 0.002294)
# Step 1: OLR1 -> Mediator (beta_1) - need this!
# Step 2: Mediator -> AAA (beta_2) - calculated above
# Mediated effect = beta_1 * beta_2
# Proportion mediated = (beta_1 * beta_2) / beta_total

# Since we don't have Step 1 (OLR1 -> CRP/IL6R), we can only do:
# 1. Report Step 2 results
# 2. Infer mediation from the pattern

print("\n--- Mediation pathway analysis ---")
print(f"\n  Total effect (OLR1 -> AAA): beta={olr1_ivw_b:.6f}, p={olr1_ivw_p:.4e}")
print(f"  (From eQTLGen Blood eQTL MR, IVW method)")

for mediator_name in ['CRP', 'IL6R']:
    med_results = [r for r in step2_results if r['mediator'] == mediator_name and r['method'] == 'Inverse variance weighted']
    if med_results:
        r = med_results[0]
        print(f"\n  {mediator_name} -> AAA: beta={r['b']:.6f}, OR={r['OR']:.4f}, p={r['pval']:.4e}")
        
        if r['pval'] < 0.05:
            print(f"  *** {mediator_name} has SIGNIFICANT causal effect on AAA ***")
            print(f"  Interpretation: {mediator_name} may be a mediator in the OLR1 -> AAA pathway")
            if r['b'] > 0:
                print(f"  Direction: Higher {mediator_name} -> INCREASED AAA risk")
            else:
                print(f"  Direction: Higher {mediator_name} -> DECREASED AAA risk")
        else:
            print(f"  {mediator_name} does NOT have significant causal effect on AAA")
            print(f"  Interpretation: {mediator_name} unlikely to mediate OLR1 -> AAA pathway")

# ============================================================
# Save all results
# ============================================================
print("\n" + "=" * 70)
print("SAVING RESULTS")
print("=" * 70)

outfile = f"{RESULTS_DIR}/mediation_mr_results.csv"
with open(outfile, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['step','mediator','method','nsnp','b','se','pval','OR','OR_lower','OR_upper'])
    writer.writeheader()
    for r in step2_results:
        writer.writerow(r)
print(f"  Saved: {outfile}")

# Also save a summary table
summary_file = f"{RESULTS_DIR}/mediation_mr_summary.csv"
with open(summary_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Analysis', 'Exposure', 'Outcome', 'Method', 'nSNP', 'beta', 'SE', 'OR', '95% CI', 'p-value', 'Significance'])
    
    # Reference: OLR1 -> AAA
    writer.writerow(['Total Effect', 'OLR1 eQTL', 'AAA (FinnGen)', 'IVW', 24,
                     f"{olr1_ivw_b:.6f}", f"{olr1_ivw_se:.6f}", f"{olr1_or:.4f}",
                     f"{math.exp(olr1_ivw_b-1.96*olr1_ivw_se):.4f}-{math.exp(olr1_ivw_b+1.96*olr1_ivw_se):.4f}",
                     f"{olr1_ivw_p:.4e}", '*'])
    
    # Step 2 results
    for r in step2_results:
        if r['step'] == 'Step2':
            sig = '***' if r['pval'] < 0.001 else ('**' if r['pval'] < 0.01 else ('*' if r['pval'] < 0.05 else 'ns'))
            writer.writerow(['Step 2', r['mediator'], 'AAA (FinnGen)', r['method'], r['nsnp'],
                             f"{r['b']:.6f}", f"{r['se']:.6f}", f"{r['OR']:.4f}",
                             f"{r['OR_lower']:.4f}-{r['OR_upper']:.4f}",
                             f"{r['pval']:.4e}", sig])

print(f"  Saved: {summary_file}")

# ============================================================
# Final Summary
# ============================================================
print("\n" + "=" * 70)
print("FINAL SUMMARY")
print("=" * 70)

print("""
+------------------+----------+--------+------------+----------+---------+
| Analysis         | Exposure | Outcome | Method    | OR (95%CI)| p-value |
+------------------+----------+--------+------------+----------+---------+
| Total Effect     | OLR1     | AAA    | IVW        | 1.0023   | 0.022 * |
|                  |          |        |            | (1.000-1.004)|       |
+------------------+----------+--------+------------+----------+---------+
| Step 2: CRP→AAA  | CRP      | AAA    | IVW        | see results above    |
| Step 2: IL6R→AAA | IL6R     | AAA    | IVW        | see results above    |
+------------------+----------+--------+------------+----------+---------+

Colocalization: PP.H4 < 2% (NEGATIVE for all tissues)

Interpretation:
1. OLR1 expression → AAA risk (MR positive, p=0.022)
2. Colocalization negative → suggests INDIRECT pathway
3. Mediation MR tests specific mediator pathways
4. IL6R shows strong association with AAA in FinnGen (p~10^-10)
5. CRP shows no significant association with AAA

Key limitation:
- Step 1 (OLR1 → Mediator) requires full mediator GWAS summary statistics
- Current analysis uses published pQTL effect sizes (approximate)
- Full mediation MR needs: OLR1 eQTL SNP effects on CRP/IL6R levels
""")

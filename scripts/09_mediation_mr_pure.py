"""
Mediation MR Analysis - Pure Python (no numpy/scipy required)
==============================================================
Step 2: Mediator (CRP/IL6R) -> AAA (FinnGen)
Step 1: OLR1 -> Mediator (GWAS Catalog lookup + literature)
"""
import os, csv, gzip, json, math, sys, time
import urllib.request, ssl
from statistics import NormalDist

OUTDIR = "D:/LOX1_EVAR_MR"
RESULTS_DIR = os.path.join(OUTDIR, "results", "tables")
FIGURES_DIR = os.path.join(OUTDIR, "results", "figures")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# ============================================================
# Statistics helpers (pure Python)
# ============================================================
nd = NormalDist()

def norm_sf(x):
    """Survival function (1 - CDF) of standard normal"""
    return nd.sf(abs(x)) * 2  # Two-tailed

def chi2_sf(x, df):
    """Survival function of chi-square distribution"""
    if df <= 0:
        return float('nan')
    # Use the regularized upper incomplete gamma function
    # P(X > x) for chi-square with df degrees of freedom
    return math.gamma_upper(df/2, x/2) if hasattr(math, 'gamma_upper') else _gamma_q(df/2, x/2)

def _gamma_q(s, x):
    """Regularized upper incomplete gamma function Q(s, x) = Gamma(s, x) / Gamma(s)"""
    # Using continued fraction expansion (Numerical Recipes)
    if x < 0 or s <= 0:
        return float('nan')
    if x == 0:
        return 1.0
    if x < s + 1:
        # Use series representation
        return 1.0 - _gamma_p_series(s, x)
    else:
        # Use continued fraction
        return _gamma_q_cf(s, x)

def _gamma_p_series(s, x):
    """Series representation of regularized lower incomplete gamma P(s, x)"""
    if x < 0:
        return 0.0
    if x == 0:
        return 0.0
    # ln(x^s * e^(-x) / Gamma(s)) = s*ln(x) - x - lgamma(s)
    ln_pref = s * math.log(x) - x - math.lgamma(s)
    sum_val = 1.0
    term = 1.0
    for n in range(1, 1000):
        term *= x / (s + n - 1)
        sum_val += term
        if abs(term) < 1e-15 * abs(sum_val):
            break
    return math.exp(ln_pref) * sum_val

def _gamma_q_cf(s, x):
    """Continued fraction representation of Q(s, x)"""
    # Lentz's method
    tiny = 1e-30
    ln_pref = s * math.log(x) - x - math.lgamma(s)
    
    b = x + 1 - s
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, 1000):
        an = -i * (i - s)
        b += 2
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    return math.exp(ln_pref) * h

# ============================================================
# MR Methods (pure Python)
# ============================================================
def wald_ratio(b_exp, se_exp, b_out, se_out):
    b = b_out / b_exp
    se = abs(b) * math.sqrt((se_exp/b_exp)**2 + (se_out/b_out)**2)
    pval = norm_sf(b/se)
    return {'method': 'Wald ratio', 'b': b, 'se': se, 'pval': pval, 'nsnp': 1}

def ivw(b_exp, se_exp, b_out, se_out):
    n = len(b_exp)
    # Weights = b_exp^2 / se_out^2
    weights = [b_exp[i]**2 / se_out[i]**2 for i in range(n)]
    
    # IVW: sum(w * b_exp * b_out) / sum(w * b_exp^2)
    num = sum(weights[i] * b_exp[i] * b_out[i] for i in range(n))
    den = sum(weights[i] * b_exp[i]**2 for i in range(n))
    b_ivw = num / den
    
    # Basic SE
    se_basic = math.sqrt(1.0 / sum(weights))
    
    # Cochran's Q
    b_ratio = [b_out[i] / b_exp[i] for i in range(n)]
    se_ratio = [se_out[i] / abs(b_exp[i]) for i in range(n)]
    Q = sum(((b_ratio[i] - b_ivw) / se_ratio[i])**2 for i in range(n))
    df = n - 1
    
    # Overdispersion
    if Q > df and df > 0:
        phi = Q / df
    else:
        phi = 1.0
    se_ivw = se_basic * math.sqrt(phi)
    
    pval = norm_sf(b_ivw / se_ivw)
    Q_pval = chi2_sf(Q, df) if df > 0 else float('nan')
    
    return {
        'method': 'Inverse variance weighted', 'b': b_ivw, 'se': se_ivw,
        'pval': pval, 'nsnp': n, 'Q': Q, 'Q_df': df, 'Q_pval': Q_pval
    }

def mr_egger(b_exp, se_exp, b_out, se_out):
    n = len(b_exp)
    # Weighted regression: b_out ~ a + b * b_exp, weights = 1/se_out^2
    w = [1.0 / se_out[i]**2 for i in range(n)]
    
    # Weighted least squares
    Sw = sum(w)
    Swx = sum(w[i] * b_exp[i] for i in range(n))
    Swy = sum(w[i] * b_out[i] for i in range(n))
    Swxx = sum(w[i] * b_exp[i]**2 for i in range(n))
    Swxy = sum(w[i] * b_exp[i] * b_out[i] for i in range(n))
    
    det = Sw * Swxx - Swx * Swx
    if abs(det) < 1e-30:
        return {'method': 'MR Egger', 'b': float('nan'), 'se': float('nan'),
                'pval': float('nan'), 'nsnp': n, 'intercept': float('nan'),
                'intercept_se': float('nan'), 'intercept_pval': float('nan')}
    
    intercept = (Swxx * Swy - Swx * Swxy) / det
    slope = (Sw * Swxy - Swx * Swy) / det
    
    # Residuals and variance
    resid = [b_out[i] - intercept - slope * b_exp[i] for i in range(n)]
    ssr = sum(w[i] * resid[i]**2 for i in range(n))
    df = n - 2
    sigma2 = ssr / df if df > 0 else 1.0
    
    se_intercept = math.sqrt(sigma2 * Swxx / det)
    se_slope = math.sqrt(sigma2 * Sw / det)
    
    pval = norm_sf(slope / se_slope)
    int_pval = norm_sf(intercept / se_intercept)
    
    return {
        'method': 'MR Egger', 'b': slope, 'se': se_slope, 'pval': pval, 'nsnp': n,
        'intercept': intercept, 'intercept_se': se_intercept, 'intercept_pval': int_pval
    }

def weighted_median(b_exp, se_exp, b_out, se_out):
    n = len(b_exp)
    b_ratio = [b_out[i] / b_exp[i] for i in range(n)]
    weights = [b_exp[i]**2 / se_out[i]**2 for i in range(n)]
    w_sum = sum(weights)
    w_norm = [w / w_sum for w in weights]
    
    # Sort
    order = sorted(range(n), key=lambda i: b_ratio[i])
    b_sorted = [b_ratio[i] for i in order]
    w_sorted = [w_norm[i] for i in order]
    
    cumsum = 0
    b_wm = b_sorted[-1]
    for i in range(n):
        cumsum += w_sorted[i]
        if cumsum >= 0.5:
            b_wm = b_sorted[i]
            break
    
    # Bootstrap SE
    import random
    random.seed(42)
    n_boot = 1000
    boot_vals = []
    for _ in range(n_boot):
        idx = [random.randint(0, n-1) for _ in range(n)]
        b_boot = [b_out[i] / b_exp[i] for i in idx]
        w_boot = [b_exp[i]**2 / se_out[i]**2 for i in idx]
        ws = sum(w_boot)
        wn = [w / ws for w in w_boot]
        order_b = sorted(range(n), key=lambda i: b_boot[i])
        b_s = [b_boot[i] for i in order_b]
        w_s = [wn[i] for i in order_b]
        cs = 0
        bv = b_s[-1]
        for i in range(n):
            cs += w_s[i]
            if cs >= 0.5:
                bv = b_s[i]
                break
        boot_vals.append(bv)
    
    se_wm = (sum((v - sum(boot_vals)/len(boot_vals))**2 for v in boot_vals) / len(boot_vals))**0.5
    pval = norm_sf(b_wm / se_wm) if se_wm > 0 else float('nan')
    
    return {'method': 'Weighted median', 'b': b_wm, 'se': se_wm, 'pval': pval, 'nsnp': n}

def run_all_mr(b_exp, se_exp, b_out, se_out, snps=None):
    results = []
    n = len(b_exp)
    if n == 0:
        return results
    if n == 1:
        results.append(wald_ratio(b_exp[0], se_exp[0], b_out[0], se_out[0]))
    else:
        results.append(ivw(b_exp, se_exp, b_out, se_out))
        results.append(mr_egger(b_exp, se_exp, b_out, se_out))
        results.append(weighted_median(b_exp, se_exp, b_out, se_out))
    return results

# ============================================================
# FinnGen extraction
# ============================================================
def complement(a):
    return {'A':'T','T':'A','G':'C','C':'G'}.get(a.upper(), a)

def extract_finngen(snp_list):
    """Extract FinnGen AAA outcomes for given SNPs"""
    path = "D:/LOX1_EVAR_MR/data/outcome/finngen_R12_I9_ABAORTANEUR.gz"
    
    targets = {}
    for s in snp_list:
        key = (int(s['chr']), int(s['pos']))
        targets.setdefault(key, []).append(s)
    
    print(f"  Searching for {len(targets)} positions in FinnGen (774MB)...")
    found = {}
    
    with gzip.open(path, 'rt') as f:
        header = f.readline().strip().split('\t')
        h = [x.lstrip('#') for x in header]
        ci = {name: i for i, name in enumerate(h)}
        
        c_chrom = ci.get('chrom', 0)
        c_pos = ci.get('pos', 1)
        c_ref = ci.get('ref', 2)
        c_alt = ci.get('alt', 3)
        c_beta = ci.get('beta', None)
        c_se = ci.get('sebeta', None)
        c_pval = ci.get('pval', None)
        c_af = ci.get('af_alt', ci.get('af', None))
        c_n = ci.get('n', None)
        
        print(f"  Columns: {h}")
        print(f"  Indices: chrom={c_chrom}, pos={c_pos}, ref={c_ref}, alt={c_alt}, beta={c_beta}, se={c_se}, pval={c_pval}, af={c_af}")
        
        cnt = 0
        for line in f:
            cnt += 1
            p = line.strip().split('\t')
            key = (int(p[c_chrom]), int(p[c_pos]))
            if key not in targets:
                continue
            
            ref = p[c_ref]
            alt = p[c_alt]
            beta = float(p[c_beta]) if c_beta is not None else float('nan')
            se = float(p[c_se]) if c_se is not None else float('nan')
            pval = float(p[c_pval]) if c_pval is not None else float('nan')
            af = float(p[c_af]) if c_af is not None else float('nan')
            n = int(float(p[c_n])) if c_n is not None else 0
            
            for s in targets[key]:
                ea = s['effect_allele']
                oa = s['other_allele']
                flip = False
                
                if ea == alt and oa == ref:
                    flip = False
                elif ea == ref and oa == alt:
                    flip = True
                elif ea == complement(alt) and oa == complement(ref):
                    flip = False
                elif ea == complement(ref) and oa == complement(alt):
                    flip = True
                else:
                    continue
                
                s['beta_outcome'] = -beta if flip else beta
                s['se_outcome'] = se
                s['pval_outcome'] = pval
                s['af_outcome'] = (1-af) if flip else af
                s['n_outcome'] = n
                found[s['SNP']] = s
                print(f"    FOUND {s['SNP']}: chr{key[0]}:{key[1]} {ref}/{alt} -> beta={s['beta_outcome']:.6f} se={s['se_outcome']:.6f} p={s['pval_outcome']:.2e}")
        
        print(f"  Scanned {cnt} lines, found {len(found)}/{len(snp_list)} SNPs")
    
    return found

# ============================================================
# GWAS Catalog lookup for mediator exposure effects
# ============================================================
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

def lookup_snp_traits(rsid, target_keywords):
    """Look up SNP associations in GWAS Catalog"""
    url = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rsid}/associations"
    try:
        req = urllib.request.Request(url, headers={'Accept': 'application/json', 'User-Agent': 'Mozilla/5.0'})
        resp = urllib.request.urlopen(req, timeout=15, context=ctx)
        data = json.loads(resp.read())
        assocs = data.get('_embedded', {}).get('associations', [])
        results = []
        for a in assocs:
            trait = ''
            loci = a.get('loci', [])
            for loc in loci:
                for rep in loc.get('authorReportedTraits', []):
                    trait = rep.get('trait', '')
                    break
            if not trait:
                trait = a.get('traitName', '')
            
            trait_lower = trait.lower()
            if any(kw in trait_lower for kw in target_keywords):
                beta = a.get('beta')
                se = a.get('standardError')
                pval = a.get('pvalue')
                odds = a.get('oddsRatio')
                acc = a.get('studyAccessionId', 'N/A')
                
                # Get risk allele
                risk_allele = ''
                if loci:
                    alleles = loci[0].get('strongestRiskAlleles', [])
                    if alleles:
                        risk_allele = alleles[0].get('riskAlleleName', '')
                
                results.append({
                    'trait': trait, 'beta': float(beta) if beta else None,
                    'se': float(se) if se else None, 'pval': float(pval) if pval else None,
                    'or': float(odds) if odds else None, 'risk_allele': risk_allele,
                    'study': acc
                })
        return results
    except:
        return []

# ============================================================
# MAIN
# ============================================================
print("=" * 70)
print("MEDIATION MR ANALYSIS")
print("=" * 70)

# --- Load OLR1 eQTL ---
print("\n[1] Loading OLR1 eQTL exposure data...")
olr1 = []
with open(f"{OUTDIR}/data/exposure/eqtlgen_OLR1_blood_eqtl.csv", 'r') as f:
    for row in csv.DictReader(f):
        olr1.append({
            'SNP': row['SNP'], 'chr': '12', 'pos': '',
            'effect_allele': row['effect_allele'], 'other_allele': row['other_allele'],
            'beta_exp': float(row['beta']),
            'se_exp': float(row['se']) if row['se'] else 0.01,
            'pval_exp': float(row['pval'])
        })
print(f"  OLR1 eQTL SNPs: {len(olr1)}")

# --- Load mediator instruments ---
print("\n[2] Loading mediator instruments...")
mediators = {}
for name in ['CRP', 'IL6R']:
    filepath = f"{OUTDIR}/data/mediators/literature_instruments_{name}.csv"
    insts = []
    with open(filepath, 'r') as f:
        for row in csv.DictReader(f):
            insts.append({
                'SNP': row['SNP'], 'chr': row['chr'], 'pos': row['pos'],
                'effect_allele': row['effect_allele'], 'other_allele': row['other_allele'],
            })
    mediators[name] = insts
    print(f"  {name}: {len(insts)} instruments")

# --- Step 2: Extract FinnGen AAA for mediator instruments ---
print("\n" + "=" * 70)
print("STEP 2: Mediator -> AAA (FinnGen)")
print("=" * 70)

all_insts = []
for name, insts in mediators.items():
    for inst in insts:
        inst['mediator'] = name
        all_insts.append(inst)

finngen_found = extract_finngen(all_insts)

# Update mediator instruments with outcomes
for name in mediators:
    matched = []
    for inst in mediators[name]:
        if inst['SNP'] in finngen_found:
            matched.append(finngen_found[inst['SNP']])
    mediators[name] = matched
    print(f"  {name}: {len(matched)} SNPs matched in FinnGen")

# --- Lookup mediator exposure effects in GWAS Catalog ---
print("\n" + "=" * 70)
print("Looking up mediator pQTL effects in GWAS Catalog...")
print("=" * 70)

# CRP keywords
crp_kw = ['c-reactive protein', 'crp']
il6_kw = ['interleukin-6', 'il-6', 'il6r', 'soluble il-6 receptor', 'interleukin 6']

mediator_exposure = {}  # SNP -> {beta_exp, se_exp, ...}

for name, insts in mediators.items():
    keywords = crp_kw if name == 'CRP' else il6_kw
    print(f"\n  {name} ({len(insts)} SNPs to look up):")
    
    for inst in insts:
        rsid = inst['SNP']
        results = lookup_snp_traits(rsid, keywords)
        if results:
            r = results[0]  # Take first matching
            inst['beta_exp'] = r['beta'] if r['beta'] else (math.log(r['or']) if r['or'] and r['or'] > 0 else None)
            inst['se_exp'] = r['se'] if r['se'] else None
            inst['pval_exp'] = r['pval']
            inst['trait_exp'] = r['trait']
            inst['study_exp'] = r['study']
            print(f"    {rsid}: {r['trait']} | beta={inst.get('beta_exp','N/A')} | p={r['pval']} | {r['study']}")
        else:
            print(f"    {rsid}: No {name} association found in GWAS Catalog")
        time.sleep(0.3)

# --- Run Step 2 MR ---
print("\n" + "=" * 70)
print("STEP 2 MR RESULTS: Mediator -> AAA")
print("=" * 70)

step2_all_results = []

for name, insts in mediators.items():
    print(f"\n--- {name} -> AAA ---")
    
    # Filter to SNPs with both exposure and outcome data
    valid = [i for i in insts if i.get('beta_exp') is not None and i.get('beta_outcome') is not None]
    
    if len(valid) == 0:
        print(f"  No SNPs with complete data. Cannot run MR.")
        # Still report FinnGen outcomes
        for inst in insts:
            if 'beta_outcome' in inst:
                print(f"  {inst['SNP']}: AAA beta={inst['beta_outcome']:.6f}, p={inst['pval_outcome']:.2e} (no exposure data)")
        continue
    
    b_exp = [i['beta_exp'] for i in valid]
    se_exp = [i['se_exp'] if i['se_exp'] else 0.01 for i in valid]
    b_out = [i['beta_outcome'] for i in valid]
    se_out = [i['se_outcome'] for i in valid]
    
    print(f"  SNPs with complete data: {len(valid)}")
    
    results = run_all_mr(b_exp, se_exp, b_out, se_out)
    
    for r in results:
        or_val = math.exp(r['b'])
        or_lower = math.exp(r['b'] - 1.96 * r['se'])
        or_upper = math.exp(r['b'] + 1.96 * r['se'])
        print(f"  {r['method']}: b={r['b']:.6f}, se={r['se']:.6f}, OR={or_val:.4f} ({or_lower:.4f}-{or_upper:.4f}), p={r['pval']:.4e}")
        
        if 'intercept' in r:
            print(f"    Egger intercept: {r['intercept']:.6f}, p={r['intercept_pval']:.4f}")
        if 'Q' in r:
            print(f"    Cochran's Q: {r['Q']:.2f} (df={r['Q_df']}), p={r['Q_pval']:.4f}")
        
        step2_all_results.append({
            'mediator': name, 'method': r['method'], 'nsnp': r['nsnp'],
            'b': r['b'], 'se': r['se'], 'pval': r['pval'],
            'OR': or_val, 'OR_lower': or_lower, 'OR_upper': or_upper,
        })

# --- Save Step 2 results ---
if step2_all_results:
    outfile = f"{RESULTS_DIR}/step2_mediator_to_aaa_mr.csv"
    with open(outfile, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['mediator','method','nsnp','b','se','pval','OR','OR_lower','OR_upper'])
        writer.writeheader()
        for r in step2_all_results:
            writer.writerow(r)
    print(f"\n  Saved: {outfile}")

# --- Save all FinnGen outcomes ---
outfile2 = f"{RESULTS_DIR}/mediator_snp_finngen_outcomes.csv"
with open(outfile2, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['mediator','SNP','chr','pos','effect_allele','other_allele',
                      'beta_outcome','se_outcome','pval_outcome','af_outcome','n_outcome',
                      'beta_exp','se_exp','pval_exp','trait_exp'])
    for name in mediators:
        for inst in mediators[name]:
            writer.writerow([
                name, inst['SNP'], inst.get('chr',''), inst.get('pos',''),
                inst.get('effect_allele',''), inst.get('other_allele',''),
                inst.get('beta_outcome',''), inst.get('se_outcome',''),
                inst.get('pval_outcome',''), inst.get('af_outcome',''),
                inst.get('n_outcome',''),
                inst.get('beta_exp',''), inst.get('se_exp',''),
                inst.get('pval_exp',''), inst.get('trait_exp','')
            ])
print(f"  Saved: {outfile2}")

# --- Summary ---
print("\n" + "=" * 70)
print("ANALYSIS SUMMARY")
print("=" * 70)

print("\n1. OLR1 -> AAA (from previous analysis):")
print("   eQTLGen IVW: OR=1.0023, p=0.022 (POSITIVE)")

print("\n2. Mediator -> AAA (Step 2):")
for name in mediators:
    insts = mediators[name]
    valid = [i for i in insts if i.get('beta_exp') is not None and i.get('beta_outcome') is not None]
    if valid:
        print(f"   {name}: {len(valid)} SNPs - see MR results above")
    else:
        print(f"   {name}: {len(insts)} FinnGen outcomes found, but no exposure betas from GWAS Catalog")
        print(f"     (Need full GWAS summary stats for proper MR)")

print("\n3. Colocalization (from previous analysis):")
print("   PP.H4 < 2% (NEGATIVE for all tissues)")

print("\n4. Interpretation:")
print("   - MR positive + coloc negative = possible indirect pathway")
print("   - Mediation MR tests: OLR1 -> [mediator] -> AAA")
print("   - Need full mediator GWAS for complete mediation analysis")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)

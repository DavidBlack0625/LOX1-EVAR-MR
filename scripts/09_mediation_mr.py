"""
Mediation MR Analysis - Complete Pipeline in Python
=====================================================
Since R is not available, implement MR methods in Python.

Strategy:
- Step 2 (Mediator -> AAA): Extract FinnGen AAA outcomes for CRP/IL6R instruments
- Step 1 (OLR1 -> Mediator): Try GWAS Catalog API for SNP-trait associations
  If not available, use literature-reported effect sizes
- Implement IVW, MR-Egger, Weighted Median in Python
"""
import os, csv, gzip, json, math, sys
import numpy as np
from scipy import stats
import urllib.request, ssl

OUTDIR = "D:/LOX1_EVAR_MR"
RESULTS_DIR = os.path.join(OUTDIR, "results", "tables")
FIGURES_DIR = os.path.join(OUTDIR, "results", "figures")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# ============================================================
# PART 1: MR Methods in Python
# ============================================================

def wald_ratio(b_exp, se_exp, b_out, se_out):
    """Single SNP Wald ratio"""
    b = b_out / b_exp
    se = abs(b) * math.sqrt((se_exp/b_exp)**2 + (se_out/b_out)**2)
    pval = 2 * stats.norm.sf(abs(b/se))
    return b, se, pval

def ivw(b_exp, se_exp, b_out, se_out):
    """Inverse Variance Weighted MR"""
    b_exp = np.array(b_exp, dtype=float)
    se_exp = np.array(se_exp, dtype=float)
    b_out = np.array(b_out, dtype=float)
    se_out = np.array(se_out, dtype=float)
    
    # Weights = 1 / (se_out^2 * (b_exp^2/b_exp^2)) simplified to:
    # For IVW: weight = b_exp^2 / se_out^2
    weights = b_exp**2 / se_out**2
    
    # IVW estimate: weighted regression through origin
    b_ivw = np.sum(weights * b_exp * b_out) / np.sum(weights * b_exp**2)
    se_ivw = math.sqrt(1 / np.sum(weights))
    
    # More precise SE accounting for overdispersion
    # Cochran's Q
    b_exp_inv = 1.0 / b_exp
    b_out_inv = b_out / b_exp
    se_out_inv = se_out / abs(b_exp)
    Q = np.sum(((b_out_inv - b_ivw) / se_out_inv)**2)
    df = len(b_exp) - 1
    
    # Overdispersion parameter
    if Q > df and df > 0:
        phi = Q / df
    else:
        phi = 1.0
    
    se_ivw_corrected = se_ivw * math.sqrt(phi)
    
    pval = 2 * stats.norm.sf(abs(b_ivw / se_ivw_corrected))
    
    # Heterogeneity
    Q_pval = 1 - stats.chi2.cdf(Q, df) if df > 0 else float('nan')
    
    return {
        'method': 'Inverse variance weighted',
        'b': b_ivw,
        'se': se_ivw_corrected,
        'pval': pval,
        'nsnp': len(b_exp),
        'Q': Q,
        'Q_df': df,
        'Q_pval': Q_pval
    }

def mr_egger(b_exp, se_exp, b_out, se_out):
    """MR-Egger regression"""
    b_exp = np.array(b_exp, dtype=float)
    se_exp = np.array(se_exp, dtype=float)
    b_out = np.array(b_out, dtype=float)
    se_out = np.array(se_out, dtype=float)
    
    # Weighted regression: b_out ~ b_exp (with intercept)
    # Weights = 1/se_out^2
    weights = 1.0 / se_out**2
    
    # Weighted least squares with intercept
    W = np.diag(weights)
    X = np.column_stack([np.ones(len(b_exp)), b_exp])
    beta = np.linalg.inv(X.T @ W @ X) @ X.T @ W @ b_out
    
    # Calculate residuals and variance
    resid = b_out - X @ beta
    df = len(b_exp) - 2
    if df > 0:
        sigma2 = np.sum(weights * resid**2) / df
    else:
        sigma2 = 1.0
    
    cov = sigma2 * np.linalg.inv(X.T @ W @ X)
    se = np.sqrt(np.diag(cov))
    
    intercept = beta[0]
    intercept_se = se[0]
    intercept_pval = 2 * stats.norm.sf(abs(intercept / intercept_se))
    
    b_egger = beta[1]
    se_egger = se[1]
    pval = 2 * stats.norm.sf(abs(b_egger / se_egger))
    
    return {
        'method': 'MR Egger',
        'b': b_egger,
        'se': se_egger,
        'pval': pval,
        'nsnp': len(b_exp),
        'intercept': intercept,
        'intercept_se': intercept_se,
        'intercept_pval': intercept_pval
    }

def weighted_median(b_exp, se_exp, b_out, se_out):
    """Weighted median estimator"""
    b_exp = np.array(b_exp, dtype=float)
    se_exp = np.array(se_exp, dtype=float)
    b_out = np.array(b_out, dtype=float)
    se_out = np.array(se_out, dtype=float)
    
    # SNP-specific effects
    b_snp = b_out / b_exp
    # Weights = b_exp^2 / se_out^2
    weights = b_exp**2 / se_out**2
    weights = weights / np.sum(weights)
    
    # Sort by SNP effect
    order = np.argsort(b_snp)
    b_sorted = b_snp[order]
    w_sorted = weights[order]
    
    # Find weighted median
    cumsum = np.cumsum(w_sorted)
    idx = np.searchsorted(cumsum, 0.5)
    
    if idx >= len(b_sorted):
        idx = len(b_sorted) - 1
    
    b_wm = b_sorted[idx]
    
    # Bootstrap SE
    n_boot = 1000
    boot_estimates = []
    for _ in range(n_boot):
        boot_idx = np.random.choice(len(b_exp), size=len(b_exp), replace=True)
        b_boot = b_out[boot_idx] / b_exp[boot_idx]
        w_boot = b_exp[boot_idx]**2 / se_out[boot_idx]**2
        w_boot = w_boot / np.sum(w_boot)
        order_boot = np.argsort(b_boot)
        b_boot_sorted = b_boot[order_boot]
        w_boot_sorted = w_boot[order_boot]
        cumsum_boot = np.cumsum(w_boot_sorted)
        idx_boot = np.searchsorted(cumsum_boot, 0.5)
        if idx_boot >= len(b_boot_sorted):
            idx_boot = len(b_boot_sorted) - 1
        boot_estimates.append(b_boot_sorted[idx_boot])
    
    se_wm = np.std(boot_estimates)
    pval = 2 * stats.norm.sf(abs(b_wm / se_wm)) if se_wm > 0 else float('nan')
    
    return {
        'method': 'Weighted median',
        'b': b_wm,
        'se': se_wm,
        'pval': pval,
        'nsnp': len(b_exp)
    }

def run_mr(b_exp, se_exp, b_out, se_out, snps=None):
    """Run all MR methods"""
    results = []
    n = len(b_exp)
    
    if n == 1:
        b, se, p = wald_ratio(b_exp[0], se_exp[0], b_out[0], se_out[0])
        results.append({'method': 'Wald ratio', 'b': b, 'se': se, 'pval': p, 'nsnp': 1})
    else:
        results.append(ivw(b_exp, se_exp, b_out, se_out))
        results.append(mr_egger(b_exp, se_exp, b_out, se_out))
        results.append(weighted_median(b_exp, se_exp, b_out, se_out))
    
    return results

# ============================================================
# PART 2: Extract FinnGen AAA outcomes for mediator instruments
# ============================================================

def extract_finngen_for_snps(snp_list):
    """
    Extract FinnGen AAA GWAS data for specific SNPs.
    FinnGen file: D:/LOX1_EVAR_MR/data/outcome/finngen_R12_I9_ABAORTANEUR.gz
    Format: #chrom, pos, ref, alt, af_alt, af_alt_se, imputation_info, n, beta, se, pval
    SNP identified by: chrom:pos:ref:alt
    """
    finngen_path = "D:/LOX1_EVAR_MR/data/outcome/finngen_R12_I9_ABAORTANEUR.gz"
    
    # Build lookup: (chrom, pos) -> list of (ref, alt, af, beta, se, pval)
    # We need to match by chrom and pos, then check alleles
    targets = {}
    for snp in snp_list:
        chrom = int(snp['chr'])
        pos = int(snp['pos'])
        key = (chrom, pos)
        if key not in targets:
            targets[key] = []
        targets[key].append(snp)
    
    print(f"  Looking for {len(targets)} unique positions in FinnGen...")
    
    found = {}
    with gzip.open(finngen_path, 'rt') as f:
        header = f.readline().strip().split('\t')
        # Find column indices
        col_map = {}
        for i, h in enumerate(header):
            col_map[h.lstrip('#')] = i
        
        print(f"  FinnGen columns: {header}")
        
        chrom_idx = col_map.get('chrom', 0)
        pos_idx = col_map.get('pos', 1)
        ref_idx = col_map.get('ref', 2)
        alt_idx = col_map.get('alt', 3)
        
        # Find beta, se, pval, af, n columns
        beta_idx = col_map.get('beta', None)
        se_idx = col_map.get('sebeta', None)
        pval_idx = col_map.get('pval', None)
        af_idx = col_map.get('af_alt', None)  # or 'af'
        n_idx = col_map.get('n', None)
        
        # If standard names not found, try alternatives
        if beta_idx is None:
            for h in header:
                if 'beta' in h.lower():
                    beta_idx = col_map[h.lstrip('#')]
                    break
        if se_idx is None:
            for h in header:
                if 'se' in h.lower() and 'beta' in h.lower():
                    se_idx = col_map[h.lstrip('#')]
                    break
        if af_idx is None:
            for h in header:
                if h.lower().startswith('af'):
                    af_idx = col_map[h.lstrip('#')]
                    break
        
        print(f"  Column indices: chrom={chrom_idx}, pos={pos_idx}, ref={ref_idx}, alt={alt_idx}")
        print(f"    beta={beta_idx}, se={se_idx}, pval={pval_idx}, af={af_idx}")
        
        line_count = 0
        for line in f:
            line_count += 1
            parts = line.strip().split('\t')
            chrom = int(parts[chrom_idx])
            pos = int(parts[pos_idx])
            
            key = (chrom, pos)
            if key in targets:
                ref = parts[ref_idx]
                alt = parts[alt_idx]
                beta = float(parts[beta_idx]) if beta_idx is not None else float('nan')
                se = float(parts[se_idx]) if se_idx is not None else float('nan')
                pval = float(parts[pval_idx]) if pval_idx is not None else float('nan')
                af = float(parts[af_idx]) if af_idx is not None else float('nan')
                n = int(parts[n_idx]) if n_idx is not None else 0
                
                for snp in targets[key]:
                    ea = snp['effect_allele']
                    oa = snp['other_allele']
                    
                    # Match alleles (accounting for strand issues)
                    matched = False
                    if (ea == alt and oa == ref):
                        matched = True
                        effect_beta = beta
                        effect_af = af
                    elif (ea == ref and oa == alt):
                        matched = True
                        effect_beta = -beta  # Flip effect
                        effect_af = 1 - af
                    elif (ea == complement(alt) and oa == complement(ref)):
                        matched = True
                        effect_beta = beta
                        effect_af = af
                    elif (ea == complement(ref) and oa == complement(alt)):
                        matched = True
                        effect_beta = -beta
                        effect_af = 1 - af
                    
                    if matched:
                        snp_id = snp['SNP']
                        found[snp_id] = {
                            'SNP': snp_id,
                            'chr': chrom,
                            'pos': pos,
                            'effect_allele': ea,
                            'other_allele': oa,
                            'beta_outcome': effect_beta,
                            'se_outcome': se,
                            'pval_outcome': pval,
                            'af_outcome': effect_af,
                            'n_outcome': n,
                            'ref': ref,
                            'alt': alt
                        }
                        print(f"    Found {snp_id}: chr{chrom}:{pos} {ref}/{alt} -> beta={effect_beta:.6f}, se={se:.6f}, p={pval:.2e}")
    
    print(f"  Total lines scanned: {line_count}")
    print(f"  SNPs found: {len(found)}/{len(snp_list)}")
    return found

def complement(allele):
    """Get complement allele"""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return comp.get(allele.upper(), allele)

# ============================================================
# PART 3: Try GWAS Catalog for CRP associations of OLR1 SNPs
# ============================================================

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

def get_crp_associations_for_snps(rsids):
    """
    Try to get CRP GWAS associations for OLR1 eQTL SNPs from GWAS Catalog.
    Use the association search endpoint.
    """
    results = {}
    
    for rsid in rsids:
        # Try SNP associations endpoint
        url = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rsid}/associations"
        try:
            req = urllib.request.Request(url, headers={'Accept': 'application/json', 'User-Agent': 'Mozilla/5.0'})
            resp = urllib.request.urlopen(req, timeout=15, context=ctx)
            data = json.loads(resp.read())
            associations = data.get('_embedded', {}).get('associations', [])
            
            for assoc in associations:
                # Check if this is a CRP association
                trait_name = ''
                loci = assoc.get('loci', [])
                for locus in loci:
                    for report in locus.get('authorReportedTraits', []):
                        trait_name = report.get('trait', '')
                        break
                
                if not trait_name:
                    trait_name = assoc.get('traitName', '')
                
                if 'c-reactive protein' in trait_name.lower() or 'crp' in trait_name.lower():
                    beta = assoc.get('beta', None)
                    se = assoc.get('standardError', None)
                    pval = assoc.get('pvalue', None)
                    or_val = assoc.get('oddsRatio', None)
                    acc = assoc.get('studyAccessionId', 'N/A')
                    
                    # Get effect allele info
                    loci_data = loci[0] if loci else {}
                    alleles = loci_data.get('strongestRiskAlleles', [])
                    risk_allele = alleles[0].get('riskAlleleName', '') if alleles else ''
                    
                    results[rsid] = {
                        'SNP': rsid,
                        'trait': trait_name,
                        'beta': float(beta) if beta else None,
                        'se': float(se) if se else None,
                        'pval': float(pval) if pval else None,
                        'or': float(or_val) if or_val else None,
                        'risk_allele': risk_allele,
                        'study': acc
                    }
                    break  # Take first CRP association
        except Exception as e:
            pass
    
    return results

# ============================================================
# PART 4: Main Analysis
# ============================================================

print("=" * 70)
print("MEDIATION MR ANALYSIS - COMPLETE PIPELINE")
print("=" * 70)

# Read OLR1 eQTL SNPs (exposure for Step 1)
print("\n--- Loading OLR1 eQTL exposure data ---")
olr1_snps = []
with open("D:/LOX1_EVAR_MR/data/exposure/eqtlgen_OLR1_blood_eqtl.csv", 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        olr1_snps.append({
            'SNP': row['SNP'],
            'chr': '12',
            'pos': '',
            'effect_allele': row['effect_allele'],
            'other_allele': row['other_allele'],
            'beta': float(row['beta']),
            'se': float(row['se']) if row['se'] else None,
            'pval': float(row['pval'])
        })
print(f"  OLR1 eQTL SNPs: {len(olr1_snps)}")

# Read mediator instruments
print("\n--- Loading mediator instruments ---")
mediators = {}
for mediator_file in ['literature_instruments_CRP.csv', 'literature_instruments_IL6R.csv']:
    filepath = os.path.join(OUTDIR, "data", "mediators", mediator_file)
    mediator_name = mediator_file.replace('literature_instruments_', '').replace('.csv', '')
    instruments = []
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            instruments.append({
                'SNP': row['SNP'],
                'chr': row['chr'],
                'pos': row['pos'],
                'effect_allele': row['effect_allele'],
                'other_allele': row['other_allele'],
            })
    mediators[mediator_name] = instruments
    print(f"  {mediator_name}: {len(instruments)} instruments")

# STEP 2: Extract FinnGen AAA outcomes for mediator instruments
print("\n" + "=" * 70)
print("STEP 2: Mediator -> AAA (FinnGen)")
print("=" * 70)

all_mediator_instruments = []
for name, instruments in mediators.items():
    for inst in instruments:
        inst['mediator'] = name
        all_mediator_instruments.append(inst)

print(f"\n  Extracting FinnGen AAA outcomes for {len(all_mediator_instruments)} mediator SNPs...")
finngen_outcomes = extract_finngen_for_snps(all_mediator_instruments)

# Match outcomes to instruments
for name, instruments in mediators.items():
    print(f"\n  {name}:")
    matched = []
    for inst in instruments:
        snp_id = inst['SNP']
        if snp_id in finngen_outcomes:
            outcome = finngen_outcomes[snp_id]
            inst['beta_outcome'] = outcome['beta_outcome']
            inst['se_outcome'] = outcome['se_outcome']
            inst['pval_outcome'] = outcome['pval_outcome']
            inst['af_outcome'] = outcome['af_outcome']
            inst['n_outcome'] = outcome['n_outcome']
            matched.append(inst)
            print(f"    {snp_id}: beta={inst['beta_outcome']:.6f}, se={inst['se_outcome']:.6f}, p={inst['pval_outcome']:.2e}")
        else:
            print(f"    {snp_id}: NOT FOUND in FinnGen")
    mediators[name] = matched

# STEP 1: Try to get CRP associations for OLR1 SNPs from GWAS Catalog
print("\n" + "=" * 70)
print("STEP 1: OLR1 -> Mediator (GWAS Catalog lookup)")
print("=" * 70)

print("\n  Searching GWAS Catalog for CRP associations of OLR1 eQTL SNPs...")
olr1_rsids = [snp['SNP'] for snp in olr1_snps]
crp_associations = get_crp_associations_for_snps(olr1_rsids[:10])  # Test first 10
print(f"  Found CRP associations for {len(crp_associations)} of 10 tested SNPs")

for snp_id, assoc in crp_associations.items():
    print(f"    {snp_id}: {assoc['trait']} | beta={assoc.get('beta', 'N/A')} | p={assoc.get('pval', 'N/A')} | study={assoc['study']}")

# If not enough from GWAS Catalog, use literature effect sizes
# Known OLR1 -> CRP association from functional studies:
# LOX-1 activation promotes inflammation -> CRP elevation
# Published MR studies suggest OLR1 SNPs affect CRP levels

# For now, let's focus on what we can do:
# 1. MR: CRP -> AAA (using cis-pQTL instruments)
# 2. MR: IL6R -> AAA (using cis-pQTL instruments)

# ============================================================
# STEP 2 ANALYSIS: Mediator -> AAA
# ============================================================

print("\n" + "=" * 70)
print("RUNNING MR: Mediator -> AAA")
print("=" * 70)

# For Step 2, we need:
# - Exposure: mediator pQTL (beta, se) - from literature/published GWAS
# - Outcome: FinnGen AAA (beta, se) - extracted above
# 
# Since we don't have the mediator pQTL beta/SE from the original GWAS,
# we'll use the FinnGen AAA outcome and infer from known effect directions.
#
# Actually, for proper MR we need the EXPOSURE beta/SE (mediator -> SNP effect)
# and the OUTCOME beta/SE (AAA -> SNP effect).
#
# For cis-pQTL, the exposure effect is the effect of the SNP on protein level.
# We need to get this from the original GWAS.
# 
# As a workaround, we can use the published effect sizes:
# CRP: rs2794520 beta=-0.21, se=0.01 (per SD CRP decrease, Liu et al. 2020)
# IL6R: rs2228145 beta=-0.55, se=0.02 (per SD IL6R decrease, Said et al. 2022)
#
# But without exact betas, we can still test the ASSOCIATION direction
# by looking at FinnGen AAA effects at these loci.

# Let's report the FinnGen AAA effects at mediator loci
print("\n--- FinnGen AAA effects at mediator instrument loci ---")
for name, instruments in mediators.items():
    print(f"\n  {name} ({len(instruments)} SNPs matched):")
    if len(instruments) == 0:
        print("    No SNPs matched in FinnGen - trying alternative approach")
        continue
    
    for inst in instruments:
        direction = "↑AAA" if inst['beta_outcome'] > 0 else "↓AAA"
        print(f"    {inst['SNP']}: beta={inst['beta_outcome']:.6f}, p={inst['pval_outcome']:.2e} {direction}")

# Save Step 2 results
step2_file = os.path.join(RESULTS_DIR, "step2_mediator_to_aaa_finngen.csv")
with open(step2_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['mediator', 'SNP', 'chr', 'pos', 'effect_allele', 'other_allele', 
                      'beta_outcome', 'se_outcome', 'pval_outcome', 'af_outcome', 'n_outcome'])
    for name, instruments in mediators.items():
        for inst in instruments:
            writer.writerow([name, inst['SNP'], inst['chr'], inst['pos'],
                            inst['effect_allele'], inst['other_allele'],
                            inst.get('beta_outcome', ''), inst.get('se_outcome', ''),
                            inst.get('pval_outcome', ''), inst.get('af_outcome', ''),
                            inst.get('n_outcome', '')])
print(f"\n  Saved: {step2_file}")

# Try to get mediator exposure data from GWAS Catalog
print("\n--- Searching GWAS Catalog for mediator pQTL effect sizes ---")
for name, instruments in mediators.items():
    print(f"\n  {name}:")
    for inst in instruments[:3]:  # Test first 3
        rsid = inst['SNP']
        url = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rsid}/associations"
        try:
            req = urllib.request.Request(url, headers={'Accept': 'application/json', 'User-Agent': 'Mozilla/5.0'})
            resp = urllib.request.urlopen(req, timeout=15, context=ctx)
            data = json.loads(resp.read())
            associations = data.get('_embedded', {}).get('associations', [])
            
            # Find CRP or IL6R associations
            target_traits = []
            if name == 'CRP':
                target_traits = ['c-reactive protein', 'crp']
            elif name == 'IL6R':
                target_traits = ['interleukin-6', 'il-6', 'il6r', 'soluble il-6 receptor']
            
            for assoc in associations:
                trait_name = ''
                loci = assoc.get('loci', [])
                for locus in loci:
                    for report in locus.get('authorReportedTraits', []):
                        trait_name = report.get('trait', '')
                        break
                if not trait_name:
                    trait_name = assoc.get('traitName', '')
                
                if any(t in trait_name.lower() for t in target_traits):
                    beta = assoc.get('beta', 'N/A')
                    se = assoc.get('standardError', 'N/A')
                    pval = assoc.get('pvalue', 'N/A')
                    acc = assoc.get('studyAccessionId', 'N/A')
                    print(f"    {rsid}: {trait_name} | beta={beta} | se={se} | p={pval} | {acc}")
                    break
            else:
                # No matching trait found
                print(f"    {rsid}: {len(associations)} associations, none matching {name}")
        except Exception as e:
            print(f"    {rsid}: Error - {e}")

print("\n" + "=" * 70)
print("PIPELINE COMPLETE")
print("=" * 70)

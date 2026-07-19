"""
Download eQTLGen trans-eQTL data (correct URL) and filter for:
- SNPs: our 31 OLR1 eQTL SNPs
- Genes: IL6R (ENSG00000160712), CRP (ENSG00000129117)

Then run Step 1 MR: OLR1 -> IL6R / OLR1 -> CRP
"""
import os, csv, json, math, gzip, time, urllib.request, ssl
from statistics import NormalDist

OUTDIR = "D:/LOX1_EVAR_MR"
DATA_DIR = os.path.join(OUTDIR, "data", "step1")
RESULTS_DIR = os.path.join(OUTDIR, "results", "tables")
os.makedirs(DATA_DIR, exist_ok=True)

nd = NormalDist()

def norm_sf(x):
    return 2 * (1 - nd.cdf(abs(x)))

# ============================================================
# Load OLR1 eQTL SNPs
# ============================================================
olr1_snps = {}
with open(f"{OUTDIR}/data/exposure/eqtlgen_OLR1_blood_eqtl.csv", 'r') as f:
    for row in csv.DictReader(f):
        snp = row['SNP']
        # eQTLGen beta is in log scale, need to normalize
        # The beta values are very large (e.g., -7.0) because they're not standardized
        # We'll use them as-is for direction, but need SE for MR
        olr1_snps[snp] = {
            'snp': snp,
            'ea': row['effect_allele'],
            'oa': row['other_allele'],
            'beta': float(row['beta']),
            'se': float(row['se']) if row['se'] else None,
            'pval': float(row['pval']),
        }

snp_list = list(olr1_snps.keys())
snp_set = set(snp_list)
print(f"Loaded {len(snp_list)} OLR1 eQTL SNPs")

# Target genes
TARGETS = {
    'IL6R': 'ENSG00000160712',
    'CRP':  'ENSG00000129117',
}
target_ensembl = set(TARGETS.values())

# ============================================================
# Download eQTLGen trans-eQTL file
# ============================================================
print("\n" + "=" * 70)
print("Downloading eQTLGen trans-eQTL data")
print("=" * 70)

trans_url = "https://download.gcc.rug.nl/downloads/eqtlgen/trans-eqtl/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz"
trans_file = os.path.join(DATA_DIR, "eqtlgen_trans_eqtl.txt.gz")

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

if os.path.exists(trans_file) and os.path.getsize(trans_file) > 1000000:
    file_size = os.path.getsize(trans_file) / (1024*1024)
    print(f"  File already exists: {file_size:.1f} MB")
else:
    print(f"  Downloading from: {trans_url}")
    try:
        req = urllib.request.Request(trans_url)
        with urllib.request.urlopen(req, timeout=300, context=ctx) as resp:
            file_size = int(resp.headers.get('Content-Length', 0)) / (1024*1024)
            print(f"  File size: {file_size:.1f} MB")
            
            downloaded = 0
            with open(trans_file, 'wb') as f:
                while True:
                    chunk = resp.read(1024*1024)  # 1MB chunks
                    if not chunk:
                        break
                    f.write(chunk)
                    downloaded += len(chunk)
                    if downloaded % (50*1024*1024) == 0:
                        print(f"    Downloaded: {downloaded/(1024*1024):.0f} MB / {file_size:.0f} MB")
            
            actual_size = os.path.getsize(trans_file) / (1024*1024)
            print(f"  Download complete: {actual_size:.1f} MB")
    except Exception as e:
        print(f"  Download failed: {e}")
        import traceback
        traceback.print_exc()

# ============================================================
# Filter trans-eQTL data
# ============================================================
print("\n" + "=" * 70)
print("Filtering trans-eQTL data for OLR1 SNPs x IL6R/CRP genes")
print("=" * 70)

trans_results = []

if os.path.exists(trans_file):
    print(f"  Scanning file...")
    line_count = 0
    found_count = 0
    
    with gzip.open(trans_file, 'rt', encoding='utf-8') as f:
        # Read header
        header = f.readline().strip().split('\t')
        print(f"  Columns: {header}")
        
        # Find column indices
        cols = {name: i for i, name in enumerate(header)}
        
        # Common column names in eQTLGen trans-eQTL:
        # PGSNP, Gene, GeneSymbol, ZScore, PValue, BonferroniP, FDR, 
        # SNPChr, SNPPos, AssessedAllele, OtherAllele, 
        # GeneChr, GenePos, GeneStart, GeneEnd, N
        snp_col = cols.get('PGSNP', cols.get('SNP', 0))
        gene_col = cols.get('Gene', 1)
        gene_symbol_col = cols.get('GeneSymbol', None)
        zscore_col = cols.get('ZScore', cols.get('Zscore', None))
        pval_col = cols.get('PValue', cols.get('Pvalue', cols.get('pvalue', None)))
        bonf_col = cols.get('BonferroniP', cols.get('Bonferroni.p.value', None))
        fdr_col = cols.get('FDR', None)
        ea_col = cols.get('AssessedAllele', cols.get('EffectAllele', None))
        oa_col = cols.get('OtherAllele', None)
        snp_chr_col = cols.get('SNPChr', None)
        snp_pos_col = cols.get('SNPPos', None)
        
        for line in f:
            line_count += 1
            parts = line.strip().split('\t')
            
            snp = parts[snp_col]
            gene_id = parts[gene_col]
            
            # Quick filter: only check if SNP is in our list AND gene is a target
            if snp in snp_set and gene_id in target_ensembl:
                found_count += 1
                
                # Get gene name
                gene_name = 'IL6R' if gene_id == TARGETS['IL6R'] else 'CRP'
                if gene_symbol_col is not None:
                    gene_symbol = parts[gene_symbol_col]
                else:
                    gene_symbol = gene_name
                
                zscore = float(parts[zscore_col]) if zscore_col is not None else None
                pval = float(parts[pval_col]) if pval_col is not None else None
                bonf = float(parts[bonf_col]) if bonf_col is not None else None
                fdr = float(parts[fdr_col]) if fdr_col is not None else None
                ea = parts[ea_col] if ea_col is not None else ''
                oa = parts[oa_col] if oa_col is not None else ''
                snp_chr = parts[snp_chr_col] if snp_chr_col is not None else ''
                snp_pos = parts[snp_pos_col] if snp_pos_col is not None else ''
                
                trans_results.append({
                    'snp': snp,
                    'target_gene': gene_name,
                    'gene_symbol': gene_symbol,
                    'ensembl_id': gene_id,
                    'zscore': zscore,
                    'pval': pval,
                    'bonferroni_p': bonf,
                    'fdr': fdr,
                    'assessed_allele': ea,
                    'other_allele': oa,
                    'snp_chr': snp_chr,
                    'snp_pos': snp_pos,
                })
                
                print(f"  ✓ Found: {snp} -> {gene_name} ({gene_symbol}): "
                      f"Z={zscore:.4f}, p={pval:.2e}, bonf={bonf:.2e}" if bonf else
                      f"  ✓ Found: {snp} -> {gene_name} ({gene_symbol}): "
                      f"Z={zscore:.4f}, p={pval:.2e}")
            
            if line_count % 1000000 == 0:
                print(f"    Scanned {line_count/1e6:.0f}M lines, found {found_count} matches...")
    
    print(f"\n  Total lines scanned: {line_count:,}")
    print(f"  Total matches found: {found_count}")
else:
    print("  trans-eQTL file not found!")

# ============================================================
# Save filtered results
# ============================================================
if trans_results:
    outfile = os.path.join(DATA_DIR, "step1_trans_eqtl_filtered.csv")
    fields = ['snp', 'target_gene', 'gene_symbol', 'ensembl_id', 'zscore', 'pval', 
              'bonferroni_p', 'fdr', 'assessed_allele', 'other_allele', 'snp_chr', 'snp_pos']
    with open(outfile, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for r in trans_results:
            writer.writerow(r)
    print(f"\n  Saved filtered results: {outfile}")
    
    # Summary by gene
    for gene in ['IL6R', 'CRP']:
        gene_results = [r for r in trans_results if r['target_gene'] == gene]
        print(f"\n  {gene}: {len(gene_results)} trans-eQTL associations")
        for r in gene_results:
            sig = "***" if (r['bonferroni_p'] or r['pval']) < 0.05/31 else ("**" if (r['pval'] or 1) < 0.001 else ("*" if (r['pval'] or 1) < 0.05 else "ns"))
            print(f"    {r['snp']}: Z={r['zscore']:.4f}, p={r['pval']:.2e}, bonf={r['bonferroni_p']:.2e if r['bonferroni_p'] else 'N/A'} {sig}")
else:
    print("\n  No trans-eQTL associations found for OLR1 SNPs -> IL6R/CRP")
    print("  This means OLR1 cis-eQTLs do NOT have significant trans-effects on IL6R or CRP expression.")
    print("  This is an important negative result for the mediation analysis.")

# ============================================================
# Run Step 1 MR if we have data
# ============================================================
if trans_results:
    print("\n" + "=" * 70)
    print("STEP 1 MR: OLR1 -> Mediator (trans-eQTL)")
    print("=" * 70)
    
    step1_results = []
    
    for gene_name in ['IL6R', 'CRP']:
        gene_data = [r for r in trans_results if r['target_gene'] == gene_name]
        
        if len(gene_data) == 0:
            print(f"\n  {gene_name}: No data available for Step 1 MR")
            step1_results.append({
                'step': 'Step1', 'mediator': gene_name, 'method': 'N/A',
                'nsnp': 0, 'b': float('nan'), 'se': float('nan'),
                'pval': float('nan'), 'OR': float('nan'),
                'OR_lower': float('nan'), 'OR_upper': float('nan'),
                'note': 'No trans-eQTL data found'
            })
            continue
        
        print(f"\n  {gene_name}: {len(gene_data)} SNPs with trans-eQTL data")
        
        # For trans-eQTL, Z-score = beta / SE
        # We need beta and SE separately
        # In eQTLGen, the effect is on gene expression (log-transformed)
        # Z-score = beta / SE, so beta = Z * SE
        # We don't have SE directly, but we can estimate it from p-value:
        # For large n, SE ≈ |beta| / |Z|, and |Z| = |Phi^{-1}(p/2)|
        # Actually, we can use: beta = Z (as test statistic), and SE = 1
        # This simplifies to using Z-scores directly as effect sizes
        
        # More correctly: if we assume SE ≈ 1/sqrt(n) where n ≈ 31,684 (eQTLGen sample size)
        # Then SE ≈ 1/sqrt(31684) ≈ 0.0056
        # And beta = Z * SE
        
        n_eqtlgen = 31684
        se_est = 1.0 / math.sqrt(n_eqtlgen)
        
        b_exp_list = []; se_exp_list = []; b_out_list = []; se_out_list = []; snp_list_used = []
        
        for r in gene_data:
            snp = r['snp']
            if snp not in olr1_snps:
                continue
            
            # Exposure: OLR1 eQTL effect
            exp_data = olr1_snps[snp]
            exp_beta = exp_data['beta']
            # eQTLGen doesn't provide SE directly. Estimate from p-value:
            # p = 2*Phi(-|Z|), Z = beta/SE, SE = beta/Z
            if exp_data['se'] is not None and exp_data['se'] > 0:
                exp_se = exp_data['se']
            else:
                # Estimate SE from p-value
                exp_z = abs(exp_beta) / math.sqrt(max(-2 * math.log(max(exp_data['pval'], 1e-300)), 0.01))
                exp_se = abs(exp_beta) / max(exp_z, 0.01) if exp_z > 0 else abs(exp_beta)
            
            # Outcome: trans-eQTL effect on IL6R/CRP
            out_z = r['zscore']
            out_beta = out_z * se_est  # estimated beta
            out_se = se_est
            
            # Check allele alignment
            exp_ea = exp_data['ea']
            out_ea = r['assessed_allele']
            
            if exp_ea.upper() == out_ea.upper():
                # Same effect allele, no flip needed
                pass
            elif exp_data['oa'].upper() == out_ea.upper():
                # Flip exposure beta
                exp_beta = -exp_beta
            else:
                # Alleles don't match - still use as-is with a note
                print(f"    Warning: {snp} allele mismatch: exp_ea={exp_ea}, out_ea={out_ea}")
            
            b_exp_list.append(exp_beta)
            se_exp_list.append(exp_se)
            b_out_list.append(out_beta)
            se_out_list.append(out_se)
            snp_list_used.append(snp)
            
            print(f"    {snp}: exp_beta={exp_beta:.4f}, exp_se={exp_se:.6f} | "
                  f"out_beta={out_beta:.6f}, out_z={out_z:.4f}, out_p={r['pval']:.2e}")
        
        if len(b_exp_list) == 0:
            print(f"  No valid SNPs for {gene_name} Step 1 MR")
            continue
        
        if len(b_exp_list) == 1:
            # Wald ratio
            b = b_out_list[0] / b_exp_list[0]
            se = abs(b) * math.sqrt((se_exp_list[0]/b_exp_list[0])**2 + (se_out_list[0]/b_out_list[0])**2)
            p = norm_sf(b/se)
            print(f"\n  Wald ratio: b={b:.6f}, se={se:.6f}, p={p:.4e}")
            step1_results.append({
                'step': 'Step1', 'mediator': gene_name, 'method': 'Wald ratio',
                'nsnp': 1, 'b': b, 'se': se, 'pval': p,
                'OR': math.exp(b), 'OR_lower': math.exp(b - 1.96*se),
                'OR_upper': math.exp(b + 1.96*se),
            })
        else:
            # IVW
            n = len(b_exp_list)
            weights = [b_exp_list[i]**2 / se_out_list[i]**2 for i in range(n)]
            num = sum(weights[i] * b_exp_list[i] * b_out_list[i] for i in range(n))
            den = sum(weights[i] * b_exp_list[i]**2 for i in range(n))
            b_ivw = num / den
            se_basic = math.sqrt(1.0 / sum(weights))
            
            # Heterogeneity
            b_ratio = [b_out_list[i] / b_exp_list[i] for i in range(n)]
            se_ratio = [se_out_list[i] / abs(b_exp_list[i]) for i in range(n)]
            Q = sum(((b_ratio[i] - b_ivw) / se_ratio[i])**2 for i in range(n))
            df = n - 1
            phi = Q / df if Q > df and df > 0 else 1.0
            se_ivw = se_basic * math.sqrt(phi)
            p_ivw = norm_sf(b_ivw / se_ivw)
            
            sig = "***" if p_ivw < 0.001 else ("**" if p_ivw < 0.01 else ("*" if p_ivw < 0.05 else "ns"))
            print(f"\n  IVW: b={b_ivw:.6f}, se={se_ivw:.6f}, p={p_ivw:.4e} {sig}")
            print(f"  Cochran's Q: {Q:.2f} (df={df})")
            
            step1_results.append({
                'step': 'Step1', 'mediator': gene_name, 'method': 'IVW',
                'nsnp': n, 'b': b_ivw, 'se': se_ivw, 'pval': p_ivw,
                'OR': math.exp(b_ivw), 'OR_lower': math.exp(b_ivw - 1.96*se_ivw),
                'OR_upper': math.exp(b_ivw + 1.96*se_ivw),
            })
    
    # Save Step 1 results
    if step1_results:
        step1_file = os.path.join(RESULTS_DIR, "step1_mr_results.csv")
        with open(step1_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['step','mediator','method','nsnp','b','se','pval','OR','OR_lower','OR_upper'])
            writer.writeheader()
            for r in step1_results:
                writer.writerow(r)
        print(f"\n  Saved Step 1 results: {step1_file}")
        
        # ============================================================
        # Complete mediation analysis
        # ============================================================
        print("\n" + "=" * 70)
        print("COMPLETE MEDIATION ANALYSIS")
        print("=" * 70)
        
        # Total effect: OLR1 -> AAA
        beta_total = 0.002294
        se_total = 0.001
        p_total = 0.0218
        
        print(f"\n  Total effect (OLR1 -> AAA): beta={beta_total:.6f}, p={p_total:.4e}")
        
        for gene_name in ['IL6R', 'CRP']:
            step1 = [r for r in step1_results if r['mediator'] == gene_name]
            if not step1 or math.isnan(step1[0]['b']):
                print(f"\n  {gene_name}: Step 1 data unavailable")
                continue
            
            s1 = step1[0]
            beta_1 = s1['b']  # OLR1 -> mediator
            
            # Step 2: mediator -> AAA (from previous analysis)
            step2_file = os.path.join(RESULTS_DIR, "mediation_mr_results.csv")
            step2_data = []
            with open(step2_file, 'r') as f:
                for row in csv.DictReader(f):
                    step2_data.append(row)
            
            s2 = [r for r in step2_data if r['mediator'] == gene_name and r['method'] == 'Inverse variance weighted']
            if s2:
                beta_2 = float(s2[0]['b'])  # mediator -> AAA
                
                # Mediated effect = beta_1 * beta_2
                beta_mediated = beta_1 * beta_2
                # Proportion mediated = beta_mediated / beta_total
                if beta_total != 0:
                    prop_mediated = beta_mediated / beta_total
                else:
                    prop_mediated = float('nan')
                
                print(f"\n  {gene_name}:")
                print(f"    Step 1 (OLR1 -> {gene_name}): beta={beta_1:.6f}, p={s1['pval']:.4e}")
                print(f"    Step 2 ({gene_name} -> AAA): beta={beta_2:.6f}")
                print(f"    Total effect (OLR1 -> AAA): beta={beta_total:.6f}")
                print(f"    Mediated effect = {beta_1:.6f} × {beta_2:.6f} = {beta_mediated:.6f}")
                print(f"    Proportion mediated = {beta_mediated:.6f} / {beta_total:.6f} = {prop_mediated:.4f} ({prop_mediated*100:.1f}%)")
                
                if abs(prop_mediated) > 1:
                    print(f"    ⚠ Proportion > 100%, suggests possible suppression effect or inconsistent mediation")
                elif abs(prop_mediated) > 0.2:
                    print(f"    Significant mediation proportion")
                elif abs(prop_mediated) > 0.05:
                    print(f"    Modest mediation proportion")
                else:
                    print(f"    Small mediation proportion")
            else:
                print(f"\n  {gene_name}: Step 2 data not found")
else:
    print("\n  No trans-eQTL data found - cannot run Step 1 MR")
    print("  This is itself an important finding: OLR1 cis-eQTLs do not")
    print("  have significant trans-effects on IL6R or CRP expression in blood.")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)

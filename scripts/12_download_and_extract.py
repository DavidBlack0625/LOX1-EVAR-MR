"""
Download FinnGen AAA GWAS with correct URL (includes /release/ path).
Then extract outcomes for mediator instrument SNPs (CRP, IL6R).
Also try to get CRP/IL6R pQTL effect sizes from GWAS Catalog.
"""
import urllib.request, ssl, os, json, csv, gzip, math, time, re

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

OUTDIR = "D:/LOX1_EVAR_MR"

# ============================================================
# 1. Download FinnGen AAA GWAS (correct URL with /release/)
# ============================================================
print("=" * 70)
print("DOWNLOADING FINNGEN R12 AAA GWAS")
print("=" * 70)

finngen_url = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_I9_ABAORTANEUR.gz"
finngen_path = os.path.join(OUTDIR, "data", "outcome", "finngen_R12_I9_ABAORTANEUR.gz")

if os.path.exists(finngen_path) and os.path.getsize(finngen_path) > 100*1024*1024:
    print(f"  Already exists: {os.path.getsize(finngen_path)/(1024*1024):.1f} MB")
else:
    print(f"  URL: {finngen_url}")
    try:
        req = urllib.request.Request(finngen_url, headers={'User-Agent': 'Mozilla/5.0'})
        resp = urllib.request.urlopen(req, timeout=30, context=ctx)
        total = int(resp.headers.get('content-length', 0))
        print(f"  File size: {total/(1024*1024):.1f} MB")
        
        downloaded = 0
        chunk = 1024 * 1024
        with open(finngen_path, 'wb') as f:
            while True:
                data = resp.read(chunk)
                if not data:
                    break
                f.write(data)
                downloaded += len(data)
                if downloaded % (100*1024*1024) < chunk:
                    pct = downloaded * 100 / total if total > 0 else 0
                    print(f"  {downloaded/(1024*1024):.0f}/{total/(1024*1024):.0f} MB ({pct:.0f}%)")
        
        actual_size = os.path.getsize(finngen_path)
        print(f"  Downloaded: {actual_size/(1024*1024):.1f} MB")
        
        if actual_size < 1000:
            print("  ERROR: File too small, download might have failed")
            with open(finngen_path, 'r') as f:
                print(f"  Content: {f.read()[:200]}")
    except Exception as e:
        print(f"  Download failed: {e}")

# ============================================================
# 2. Extract FinnGen AAA outcomes for mediator SNPs
# ============================================================
print("\n" + "=" * 70)
print("EXTRACTING FINNGEN AAA OUTCOMES FOR MEDIATOR SNPs")
print("=" * 70)

# Use GRCh38 positions from GWAS Catalog
mediator_snps = {
    'CRP': [
        {'SNP': 'rs2794520', 'chr': 1, 'pos': 159709026, 'ea': 'C', 'oa': 'T'},
        {'SNP': 'rs3093077', 'chr': 1, 'pos': 159690029, 'ea': 'C', 'oa': 'T'},
        {'SNP': 'rs1800947', 'chr': 1, 'pos': 159673627, 'ea': 'G', 'oa': 'C'},
        {'SNP': 'rs1205',    'chr': 1, 'pos': 159712443, 'ea': 'C', 'oa': 'T'},
        {'SNP': 'rs1130864', 'chr': 1, 'pos': 159713301, 'ea': 'C', 'oa': 'T'},
        {'SNP': 'rs3091244', 'chr': 1, 'pos': 159714875, 'ea': 'C', 'oa': 'T'},
        {'SNP': 'rs7552218', 'chr': 1, 'pos': 159683715, 'ea': 'A', 'oa': 'G'},
        {'SNP': 'rs12133604','chr': 1, 'pos': 159667065, 'ea': 'C', 'oa': 'T'},
    ],
    'IL6R': [
        {'SNP': 'rs2228145', 'chr': 1, 'pos': 154454494, 'ea': 'A', 'oa': 'C'},
        {'SNP': 'rs4537545', 'chr': 1, 'pos': 154446403, 'ea': 'C', 'oa': 'T'},
        {'SNP': 'rs7549338', 'chr': 1, 'pos': 154467421, 'ea': 'G', 'oa': 'A'},
        {'SNP': 'rs4129267', 'chr': 1, 'pos': 154453788, 'ea': 'C', 'oa': 'T'},
        {'SNP': 'rs2229238', 'chr': 1, 'pos': 154465420, 'ea': 'C', 'oa': 'T'},
    ]
}

# Also extract OLR1 eQTL SNPs (for Step 1 lookup)
olr1_snps = []
with open(os.path.join(OUTDIR, "data", "exposure", "eqtlgen_OLR1_blood_eqtl.csv"), 'r') as f:
    for row in csv.DictReader(f):
        olr1_snps.append(row['SNP'])

# Build target positions
all_targets = {}
for name, insts in mediator_snps.items():
    for inst in insts:
        key = (inst['chr'], inst['pos'])
        all_targets[key] = inst

print(f"  Looking for {len(all_targets)} mediator SNP positions in FinnGen...")
print(f"  (Also scanning {len(olr1_snps)} OLR1 eQTL SNPs by rsID)")

# Also need to find OLR1 SNPs - they're on chr12, positions unknown
# Let me scan the FinnGen file and match by position for mediators
# and by rsID for OLR1 (FinnGen doesn't have rsID column, only chrom:pos:ref:alt)

# First, let's check FinnGen column format
print("\n  Scanning FinnGen file...")

found = {}
olr1_region_found = []

with gzip.open(finngen_path, 'rt') as f:
    header = f.readline().strip().split('\t')
    h = [x.lstrip('#') for x in header]
    ci = {name: i for i, name in enumerate(h)}
    
    print(f"  Columns: {h}")
    
    c_chrom = ci.get('chrom', 0)
    c_pos = ci.get('pos', 1)
    c_ref = ci.get('ref', 2)
    c_alt = ci.get('alt', 3)
    c_beta = ci.get('beta', None)
    c_se = ci.get('sebeta', None)
    c_pval = ci.get('pval', None)
    c_af = ci.get('af_alt', ci.get('af', None))
    c_n = ci.get('n', None)
    
    line_count = 0
    for line in f:
        line_count += 1
        p = line.strip().split('\t')
        chrom = int(p[c_chrom])
        pos = int(p[c_pos])
        
        key = (chrom, pos)
        if key in all_targets:
            ref = p[c_ref]
            alt = p[c_alt]
            beta = float(p[c_beta]) if c_beta is not None else float('nan')
            se = float(p[c_se]) if c_se is not None else float('nan')
            pval = float(p[c_pval]) if c_pval is not None else float('nan')
            af = float(p[c_af]) if c_af is not None else float('nan')
            n = int(float(p[c_n])) if c_n is not None else 0
            
            inst = all_targets[key]
            ea = inst['ea']
            oa = inst['oa']
            
            comp = {'A':'T','T':'A','G':'C','C':'G'}
            flip = False
            
            if ea == alt and oa == ref:
                flip = False
            elif ea == ref and oa == alt:
                flip = True
            elif ea == comp.get(alt,'') and oa == comp.get(ref,''):
                flip = False
            elif ea == comp.get(ref,'') and oa == comp.get(alt,''):
                flip = True
            else:
                # Try matching with just alleles
                continue
            
            inst['beta_aaa'] = -beta if flip else beta
            inst['se_aaa'] = se
            inst['pval_aaa'] = pval
            inst['af_aaa'] = (1-af) if flip else af
            inst['n_aaa'] = n
            inst['ref_aaa'] = ref
            inst['alt_aaa'] = alt
            found[inst['SNP']] = inst
            
            direction = "risk" if inst['beta_aaa'] > 0 else "protective"
            print(f"    FOUND {inst['SNP']}: chr{chrom}:{pos} {ref}/{alt} -> beta={inst['beta_aaa']:.6f}, se={inst['se_aaa']:.6f}, p={inst['pval_aaa']:.2e} ({direction})")
        
        if line_count % 5000000 == 0:
            print(f"    ... scanned {line_count/1e6:.0f}M lines")

print(f"\n  Total lines: {line_count}")
print(f"  Mediator SNPs found: {len(found)}/{len(all_targets)}")

# ============================================================
# 3. Try GWAS Catalog for CRP/IL6R pQTL effects
# ============================================================
print("\n" + "=" * 70)
print("GWAS CATALOG: SEARCH FOR CRP/IL6R pQTL EFFECTS")
print("=" * 70)

def lookup_all_traits(rsid, keywords):
    """Get ALL associations for a SNP and filter by keywords"""
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
            
            if any(kw in trait.lower() for kw in keywords):
                beta = a.get('beta')
                se = a.get('standardError')
                pval = a.get('pvalue')
                odds = a.get('oddsRatio')
                acc = a.get('studyAccessionId', 'N/A')
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

crp_kw = ['c-reactive protein', 'crp']
il6_kw = ['interleukin-6', 'il-6', 'il6', 'soluble il-6', 'interleukin 6']

# Search for CRP and IL6R pQTL effects
for name, insts in mediator_snps.items():
    keywords = crp_kw if name == 'CRP' else il6_kw
    print(f"\n  {name} pQTL search:")
    
    for inst in insts:
        if inst['SNP'] not in found:
            continue  # Skip if not found in FinnGen
        
        results = lookup_all_traits(inst['SNP'], keywords)
        if results:
            best = max(results, key=lambda x: x['pval'] or 1)
            inst['beta_med'] = best['beta'] if best['beta'] else (math.log(best['or']) if best['or'] and best['or'] > 0 else None)
            inst['se_med'] = best['se']
            inst['pval_med'] = best['pval']
            inst['trait_med'] = best['trait']
            inst['study_med'] = best['study']
            
            b_str = f"beta={best['beta']:.4f}" if best['beta'] else (f"OR={best['or']:.4f}" if best['or'] else "N/A")
            print(f"    {inst['SNP']}: {best['trait'][:40]} | {b_str} | p={best['pval']:.2e} | {best['study']}")
        else:
            print(f"    {inst['SNP']}: No {name} association in GWAS Catalog")
        time.sleep(0.3)

# ============================================================
# 4. Save results
# ============================================================
print("\n" + "=" * 70)
print("SAVING RESULTS")
print("=" * 70)

outfile = os.path.join(OUTDIR, "results", "tables", "mediator_finngen_outcomes.csv")
with open(outfile, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['mediator', 'SNP', 'chr', 'pos', 'ea', 'oa',
                      'beta_aaa', 'se_aaa', 'pval_aaa', 'af_aaa', 'n_aaa',
                      'ref_aaa', 'alt_aaa',
                      'beta_med', 'se_med', 'pval_med', 'trait_med', 'study_med'])
    for name, insts in mediator_snps.items():
        for inst in insts:
            if inst['SNP'] in found:
                writer.writerow([
                    name, inst['SNP'], inst['chr'], inst['pos'], inst['ea'], inst['oa'],
                    inst.get('beta_aaa', ''), inst.get('se_aaa', ''), inst.get('pval_aaa', ''),
                    inst.get('af_aaa', ''), inst.get('n_aaa', ''),
                    inst.get('ref_aaa', ''), inst.get('alt_aaa', ''),
                    inst.get('beta_med', ''), inst.get('se_med', ''),
                    inst.get('pval_med', ''), inst.get('trait_med', ''), inst.get('study_med', '')
                ])
print(f"  Saved: {outfile}")

# Summary
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
for name, insts in mediator_snps.items():
    aaa = sum(1 for i in insts if 'beta_aaa' in i)
    med = sum(1 for i in insts if 'beta_med' in i)
    both = sum(1 for i in insts if 'beta_aaa' in i and 'beta_med' in i)
    print(f"  {name}: {len(insts)} SNPs | AAA found: {aaa} | Mediator found: {med} | Both: {both}")

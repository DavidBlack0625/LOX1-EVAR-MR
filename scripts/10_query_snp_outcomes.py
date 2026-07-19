"""
Query FinnGen API and GWAS Catalog for AAA outcomes of mediator instrument SNPs.
Also look up mediator pQTL effects from GWAS Catalog.
"""
import os, csv, json, math, time, sys
import urllib.request, ssl

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

OUTDIR = "D:/LOX1_EVAR_MR"
RESULTS_DIR = os.path.join(OUTDIR, "results", "tables")

def api_get(url, timeout=30):
    req = urllib.request.Request(url, headers={
        'Accept': 'application/json',
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'
    })
    resp = urllib.request.urlopen(req, timeout=timeout, context=ctx)
    return json.loads(resp.read())

# ============================================================
# Mediator instruments
# ============================================================
instruments = {
    'CRP': [
        {'SNP': 'rs2794520', 'chr': 1, 'pos': 159681684, 'ea': 'C', 'oa': 'T'},
        {'SNP': 'rs3093077', 'chr': 1, 'pos': 159690029, 'ea': 'C', 'oa': 'T'},
        {'SNP': 'rs1800947', 'chr': 1, 'pos': 159673627, 'ea': 'G', 'oa': 'C'},
        {'SNP': 'rs1205',    'chr': 1, 'pos': 159678990, 'ea': 'C', 'oa': 'T'},
        {'SNP': 'rs1130864', 'chr': 1, 'pos': 159680682, 'ea': 'C', 'oa': 'T'},
        {'SNP': 'rs3091244', 'chr': 1, 'pos': 159671067, 'ea': 'C', 'oa': 'T'},
        {'SNP': 'rs7552218', 'chr': 1, 'pos': 159683715, 'ea': 'A', 'oa': 'G'},
        {'SNP': 'rs12133604','chr': 1, 'pos': 159667065, 'ea': 'C', 'oa': 'T'},
    ],
    'IL6R': [
        {'SNP': 'rs2228145', 'chr': 1, 'pos': 154456783, 'ea': 'A', 'oa': 'C'},
        {'SNP': 'rs4537545', 'chr': 1, 'pos': 154461268, 'ea': 'C', 'oa': 'T'},
        {'SNP': 'rs7549338', 'chr': 1, 'pos': 154467421, 'ea': 'G', 'oa': 'A'},
        {'SNP': 'rs4129267', 'chr': 1, 'pos': 154426871, 'ea': 'C', 'oa': 'T'},
        {'SNP': 'rs2229238', 'chr': 1, 'pos': 154429673, 'ea': 'C', 'oa': 'T'},
    ]
}

# ============================================================
# 1. FinnGen API - Query individual variants for AAA
# ============================================================
print("=" * 70)
print("1. QUERYING FINNGEN R12 API FOR AAA OUTCOMES")
print("=" * 70)

phenocode = "I9_ABAORTANEUR"

for name, insts in instruments.items():
    print(f"\n  {name} instruments:")
    for inst in insts:
        chrom = inst['chr']
        pos = inst['pos']
        
        # Try different API URL formats
        urls = [
            f"https://r12.finngen.fi/api/individual_variant/{phenocode}/{chrom}_{pos}",
            f"https://r12.finngen.fi/api/individual_variant/{phenocode}/{chrom}:{pos}",
        ]
        
        found = False
        for url in urls:
            try:
                data = api_get(url, timeout=15)
                if data and isinstance(data, dict):
                    # Parse the result
                    beta = data.get('beta', None)
                    se = data.get('sebeta', data.get('se', None))
                    pval = data.get('pval', None)
                    af = data.get('af_alt', data.get('af', None))
                    n = data.get('n', None)
                    ref = data.get('ref', '')
                    alt = data.get('alt', '')
                    
                    if beta is not None:
                        # Determine if we need to flip
                        ea = inst['ea']
                        oa = inst['oa']
                        flip = False
                        comp = {'A':'T','T':'A','G':'C','C':'G'}
                        
                        if ea == alt and oa == ref:
                            flip = False
                        elif ea == ref and oa == alt:
                            flip = True
                        elif ea == comp.get(alt,'') and oa == comp.get(ref,''):
                            flip = False
                        elif ea == comp.get(ref,'') and oa == comp.get(alt,''):
                            flip = True
                        else:
                            print(f"    {inst['SNP']} (chr{chrom}:{pos}): allele mismatch ea={ea}/oa={oa} vs ref={ref}/alt={alt}")
                            continue
                        
                        inst['beta_aaa'] = -float(beta) if flip else float(beta)
                        inst['se_aaa'] = float(se) if se else None
                        inst['pval_aaa'] = float(pval) if pval else None
                        inst['af_aaa'] = (1-float(af)) if flip and af else (float(af) if af else None)
                        inst['n_aaa'] = int(n) if n else None
                        inst['ref_aaa'] = ref
                        inst['alt_aaa'] = alt
                        
                        direction = "risk" if inst['beta_aaa'] > 0 else "protective"
                        print(f"    {inst['SNP']}: chr{chrom}:{pos} {ref}/{alt} -> beta={inst['beta_aaa']:.6f}, se={inst['se_aaa']:.6f}, p={inst['pval_aaa']:.2e} ({direction})")
                        found = True
                        break
            except Exception as e:
                pass
        
        if not found:
            print(f"    {inst['SNP']} (chr{chrom}:{pos}): NOT FOUND via API")
        
        time.sleep(0.5)

# ============================================================
# 2. GWAS Catalog - Look up mediator pQTL effects
# ============================================================
print("\n" + "=" * 70)
print("2. GWAS CATALOG: MEDIATOR pQTL EFFECTS")
print("=" * 70)

def lookup_gwas_catalog(rsid, target_keywords):
    """Look up SNP associations in GWAS Catalog"""
    url = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rsid}/associations"
    try:
        data = api_get(url, timeout=15)
        assocs = data.get('_embedded', {}).get('associations', [])
        results = []
        for a in assocs:
            trait = ''
            loci = a.get('loci', [])
            for loc in loci:
                for rep in loc.get('authorReportedTraits', []):
                    trait = rep.get('trait', '')
                    break
            
            trait_lower = trait.lower()
            if any(kw in trait_lower for kw in target_keywords):
                beta = a.get('beta')
                se = a.get('standardError')
                pval = a.get('pvalue')
                odds = a.get('oddsRatio')
                acc = a.get('studyAccessionId', 'N/A')
                
                risk_allele = ''
                risk_allele_freq = None
                if loci:
                    alleles = loci[0].get('strongestRiskAlleles', [])
                    if alleles:
                        risk_allele = alleles[0].get('riskAlleleName', '')
                        risk_allele_freq = alleles[0].get('riskAlleleFrequency', None)
                
                results.append({
                    'trait': trait, 'beta': float(beta) if beta else None,
                    'se': float(se) if se else None, 'pval': float(pval) if pval else None,
                    'or': float(odds) if odds else None, 'risk_allele': risk_allele,
                    'risk_allele_freq': float(risk_allele_freq) if risk_allele_freq else None,
                    'study': acc
                })
        return results
    except:
        return []

crp_kw = ['c-reactive protein', 'crp']
il6_kw = ['interleukin-6', 'il-6', 'il6r', 'soluble il-6', 'interleukin 6']

for name, insts in instruments.items():
    keywords = crp_kw if name == 'CRP' else il6_kw
    print(f"\n  {name} pQTL lookup:")
    
    for inst in insts:
        results = lookup_gwas_catalog(inst['SNP'], keywords)
        if results:
            # Take the most significant association
            best = max(results, key=lambda x: x['pval'] or 1)
            inst['beta_mediator'] = best['beta'] if best['beta'] else (math.log(best['or']) if best['or'] and best['or'] > 0 else None)
            inst['se_mediator'] = best['se']
            inst['pval_mediator'] = best['pval']
            inst['trait_mediator'] = best['trait']
            inst['study_mediator'] = best['study']
            inst['risk_allele_mediator'] = best['risk_allele']
            
            beta_str = f"beta={best['beta']:.4f}" if best['beta'] else f"OR={best['or']:.4f}"
            print(f"    {inst['SNP']}: {best['trait'][:40]} | {beta_str} | p={best['pval']:.2e} | {best['study']}")
        else:
            print(f"    {inst['SNP']}: No {name} association in GWAS Catalog")
        
        time.sleep(0.3)

# ============================================================
# 3. Also look up AAA associations in GWAS Catalog
# ============================================================
print("\n" + "=" * 70)
print("3. GWAS CATALOG: AAA ASSOCIATIONS FOR MEDIATOR SNPs")
print("=" * 70)

aaa_kw = ['aortic aneurysm', 'abdominal aortic', 'aaa', 'aneurysm']

for name, insts in instruments.items():
    print(f"\n  {name} SNPs -> AAA:")
    for inst in insts:
        results = lookup_gwas_catalog(inst['SNP'], aaa_kw)
        if results:
            for r in results[:2]:
                beta_str = f"beta={r['beta']:.4f}" if r['beta'] else (f"OR={r['or']:.4f}" if r['or'] else "N/A")
                print(f"    {inst['SNP']}: {r['trait'][:50]} | {beta_str} | p={r['pval']:.2e} | {r['study']}")
        else:
            print(f"    {inst['SNP']}: No AAA association in GWAS Catalog")
        time.sleep(0.3)

# ============================================================
# 4. Save all results
# ============================================================
print("\n" + "=" * 70)
print("4. SAVING RESULTS")
print("=" * 70)

outfile = os.path.join(RESULTS_DIR, "mediator_snp_lookup_results.csv")
with open(outfile, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['mediator', 'SNP', 'chr', 'pos', 'effect_allele', 'other_allele',
                      'beta_aaa_finngen', 'se_aaa_finngen', 'pval_aaa_finngen', 'af_aaa_finngen',
                      'ref_aaa', 'alt_aaa',
                      'beta_mediator', 'se_mediator', 'pval_mediator', 'trait_mediator', 'study_mediator',
                      'risk_allele_mediator'])
    for name, insts in instruments.items():
        for inst in insts:
            writer.writerow([
                name, inst['SNP'], inst['chr'], inst['pos'], inst['ea'], inst['oa'],
                inst.get('beta_aaa', ''), inst.get('se_aaa', ''), inst.get('pval_aaa', ''),
                inst.get('af_aaa', ''), inst.get('ref_aaa', ''), inst.get('alt_aaa', ''),
                inst.get('beta_mediator', ''), inst.get('se_mediator', ''),
                inst.get('pval_mediator', ''), inst.get('trait_mediator', ''),
                inst.get('study_mediator', ''), inst.get('risk_allele_mediator', '')
            ])
print(f"  Saved: {outfile}")

# Summary
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

for name, insts in instruments.items():
    aaa_found = sum(1 for i in insts if 'beta_aaa' in i)
    med_found = sum(1 for i in insts if 'beta_mediator' in i)
    both = sum(1 for i in insts if 'beta_aaa' in i and 'beta_mediator' in i)
    print(f"  {name}: {len(insts)} total | AAA outcomes: {aaa_found} | Mediator effects: {med_found} | Both: {both}")

print("\nDone!")

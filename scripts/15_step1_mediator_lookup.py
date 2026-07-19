"""
Step 1 Mediation MR: OLR1 -> IL6R/CRP
======================================
Challenge: OLR1 cis-eQTLs are on chr12, IL6R/CRP are on chr1.
Need trans-eQTL or genome-wide pQTL data.

Approach 1: eQTLGen trans-eQTL download + filter
Approach 2: GWAS Catalog API query by SNP
Approach 3: IEU OpenGWAS API for IL6R pQTL
"""
import os, csv, json, math, gzip, time, urllib.request, urllib.error, ssl

OUTDIR = "D:/LOX1_EVAR_MR"
DATA_DIR = os.path.join(OUTDIR, "data", "step1")
os.makedirs(DATA_DIR, exist_ok=True)
RESULTS_DIR = os.path.join(OUTDIR, "results", "tables")

# Load OLR1 eQTL SNPs
olr1_snps = {}
with open(f"{OUTDIR}/data/exposure/eqtlgen_OLR1_blood_eqtl.csv", 'r') as f:
    for row in csv.DictReader(f):
        snp = row['SNP']
        olr1_snps[snp] = {
            'snp': snp,
            'ea': row['effect_allele'],
            'oa': row['other_allele'],
            'beta': float(row['beta']),
            'se': float(row['se']) if row['se'] else 0.3,  # eQTLGen SE not provided, estimate
            'pval': float(row['pval']),
        }

snp_list = list(olr1_snps.keys())
print(f"Loaded {len(snp_list)} OLR1 eQTL SNPs from eQTLGen")
print(f"First 5: {snp_list[:5]}")

# Target genes
TARGETS = {
    'IL6R': {'ensembl': 'ENSG00000160712', 'chr': 1},
    'CRP':  {'ensembl': 'ENSG00000129117', 'chr': 1},
}

# ============================================================
# Approach 1: eQTLGen trans-eQTL download
# ============================================================
print("\n" + "=" * 70)
print("APPROACH 1: eQTLGen trans-eQTL data")
print("=" * 70)

trans_url = "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/trans-eqtl/trans-eQTL_significant_20181029.txt.gz"
trans_file = os.path.join(DATA_DIR, "eqtlgen_trans_eqtl.txt.gz")

# Check if file already exists
if os.path.exists(trans_file):
    file_size = os.path.getsize(trans_file) / (1024*1024)
    print(f"  File already exists: {file_size:.1f} MB")
else:
    # Try to download - first check size with HEAD request
    try:
        req = urllib.request.Request(trans_url, method='HEAD')
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        with urllib.request.urlopen(req, timeout=30, context=ctx) as resp:
            file_size = int(resp.headers.get('Content-Length', 0)) / (1024*1024)
            print(f"  File size: {file_size:.1f} MB")
            
            if file_size < 500:  # Less than 500MB, download
                print(f"  Downloading trans-eQTL data...")
                urllib.request.urlretrieve(trans_url, trans_file)
                actual_size = os.path.getsize(trans_file) / (1024*1024)
                print(f"  Downloaded: {actual_size:.1f} MB")
            else:
                print(f"  File too large ({file_size:.0f} MB), will try streaming approach")
    except Exception as e:
        print(f"  HEAD request failed: {e}")
        # Try direct download anyway
        print(f"  Attempting direct download...")
        try:
            ctx = ssl.create_default_context()
            ctx.check_hostname = False
            ctx.verify_mode = ssl.CERT_NONE
            with urllib.request.urlopen(trans_url, timeout=60, context=ctx) as resp:
                with open(trans_file, 'wb') as f:
                    while True:
                        chunk = resp.read(1024*1024)
                        if not chunk:
                            break
                        f.write(chunk)
                actual_size = os.path.getsize(trans_file) / (1024*1024)
                print(f"  Downloaded: {actual_size:.1f} MB")
        except Exception as e2:
            print(f"  Download failed: {e2}")

# If file exists, filter for our SNPs and target genes
trans_results = {}
if os.path.exists(trans_file):
    print(f"\n  Filtering trans-eQTL data for OLR1 SNPs x IL6R/CRP genes...")
    target_genes = set(TARGETS[t]['ensembl'] for t in TARGETS)
    snp_set = set(snp_list)
    
    found = 0
    with gzip.open(trans_file, 'rt', encoding='utf-8') as f:
        header = f.readline().strip().split('\t')
        print(f"  Columns: {header}")
        
        # Find column indices
        cols = {name: i for i, name in enumerate(header)}
        snp_col = cols.get('SNP', 0)
        gene_col = cols.get('Gene', cols.get('GeneSymbol', 1))
        beta_col = cols.get('Zscore', cols.get('beta', cols.get('beta_meta', None)))
        pval_col = cols.get('Pvalue', cols.get('pval', None))
        
        for line in f:
            parts = line.strip().split('\t')
            snp = parts[snp_col]
            gene = parts[gene_col]
            
            if snp in snp_set and gene in target_genes:
                found += 1
                gene_name = 'IL6R' if gene == 'ENSG00000160712' else 'CRP'
                zscore = float(parts[beta_col]) if beta_col is not None else None
                pval = float(parts[pval_col]) if pval_col is not None else None
                
                # Convert Z-score to beta and SE
                if zscore is not None and pval is not None:
                    from statistics import NormalDist
                    nd = NormalDist()
                    # For trans-eQTL, beta = Z * SE, but we don't have SE directly
                    # Use |Z| as the test statistic
                    trans_results[(snp, gene_name)] = {
                        'snp': snp, 'gene': gene_name, 'ensembl': gene,
                        'zscore': zscore, 'pval': pval,
                    }
                    print(f"    Found: {snp} -> {gene_name}: Z={zscore:.4f}, p={pval:.2e}")
    
    print(f"\n  Total trans-eQTL matches: {found}")
    if found > 0:
        print(f"  Unique SNP-gene pairs: {len(trans_results)}")
else:
    print("  trans-eQTL file not available")

# ============================================================
# Approach 2: GWAS Catalog API - query each SNP
# ============================================================
print("\n" + "=" * 70)
print("APPROACH 2: GWAS Catalog API (query SNPs for IL6R/CRP associations)")
print("=" * 70)

gwas_catalog_results = []
api_base = "https://www.ebi.ac.uk/gwas/rest/api"

# Also search for IL6R/CRP studies in GWAS Catalog
# First, find studies with IL6R or CRP as trait
for trait_term, gene_name in [("interleukin+6+receptor", "IL6R"), ("C-reactive+protein", "CRP")]:
    search_url = f"{api_base}/studies/search/findByDiseaseTrait?diseaseTrait={trait_term}"
    try:
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        req = urllib.request.Request(search_url, headers={'Accept': 'application/json'})
        with urllib.request.urlopen(req, timeout=30, context=ctx) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            studies = data.get('_embedded', {}).get('studies', [])
            print(f"\n  Found {len(studies)} studies for '{trait_term}'")
            for s in studies[:5]:
                acc = s.get('accessionId', 'N/A')
                trait = s.get('diseaseTrait', {}).get('trait', 'N/A')
                sample_size = s.get('initialSampleSize', 'N/A')
                print(f"    {acc}: {trait} ({sample_size})")
    except Exception as e:
        print(f"  Study search failed for '{trait_term}': {e}")

# Query each SNP for associations
print(f"\n  Querying {len(snp_list)} SNPs via GWAS Catalog API...")
snp_query_results = []

for i, snp in enumerate(snp_list):
    if (i+1) % 5 == 0:
        print(f"    Progress: {i+1}/{len(snp_list)}...")
    
    snp_url = f"{api_base}/singleNucleotidePolymorphisms/{snp}/associations"
    try:
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        req = urllib.request.Request(snp_url, headers={'Accept': 'application/json'})
        with urllib.request.urlopen(req, timeout=20, context=ctx) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            associations = data.get('_embedded', {}).get('associations', [])
            
            for assoc in associations:
                trait_names = assoc.get('efoTraits', [])
                trait_labels = [t.get('trait', '') for t in trait_names]
                trait_str = '; '.join(trait_labels).lower()
                
                # Also check the study trait
                study = assoc.get('study', {})
                study_trait = study.get('diseaseTrait', {}).get('trait', '').lower() if study else ''
                
                # Check if this association is related to IL6R or CRP
                all_trait_text = trait_str + ' ' + study_trait
                
                if any(kw in all_trait_text for kw in ['interleukin-6 receptor', 'il-6 receptor', 'il6r', 'soluble il-6', 'sil6r']):
                    gene_target = 'IL6R'
                elif any(kw in all_trait_text for kw in ['c-reactive protein', 'crp', 'hs-crp', 'hscrp']):
                    gene_target = 'CRP'
                else:
                    continue
                
                # Extract effect size
                pval = assoc.get('pvalue', None)
                pval_desc = assoc.get('pvalueDescription', '')
                or_val = assoc.get('orPerCopyNum', None)
                beta_val = assoc.get('betaNum', None)
                beta_unit = assoc.get('betaUnit', '')
                se_val = assoc.get('standardError', None)
                
                # Risk allele
                risk_allele = assoc.get('riskAlleleName', '')
                loci = assoc.get('loci', [])
                if loci:
                    strongest_snp = loci[0].get('strongestRiskAlleles', [{}])[0] if loci[0].get('strongestRiskAlleles') else {}
                    risk_allele = strongest_snp.get('riskAlleleName', risk_allele)
                
                snp_query_results.append({
                    'snp': snp,
                    'target': gene_target,
                    'trait': study_trait,
                    'pval': pval,
                    'or': or_val,
                    'beta': beta_val,
                    'beta_unit': beta_unit,
                    'se': se_val,
                    'risk_allele': risk_allele,
                    'study_accession': study.get('accessionId', '') if study else '',
                })
                
                print(f"    {snp} -> {gene_target}: p={pval}, beta={beta_val} {beta_unit}, OR={or_val}, trait='{study_trait}'")
        
        time.sleep(0.3)  # Rate limiting
        
    except urllib.error.HTTPError as e:
        if e.code == 404:
            pass  # SNP not in GWAS Catalog, skip
        else:
            print(f"    {snp}: HTTP {e.code}")
    except Exception as e:
        print(f"    {snp}: Error - {e}")

print(f"\n  GWAS Catalog results: {len(snp_query_results)} associations found")

# ============================================================
# Approach 3: Try IEU OpenGWAS API
# ============================================================
print("\n" + "=" * 70)
print("APPROACH 3: IEU OpenGWAS API")
print("=" * 70)

# Search for IL6R and CRP GWAS
for trait in ["IL6R", "interleukin-6 receptor", "CRP", "C-reactive protein"]:
    search_url = f"https://gwas-api.mrcieu.ac.uk/gwasinfo?trait={trait}"
    try:
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        req = urllib.request.Request(search_url, headers={'Accept': 'application/json'})
        with urllib.request.urlopen(req, timeout=15, context=ctx) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            if data:
                for id_, info in data.items():
                    print(f"  Found: {id_} - {info.get('trait', 'N/A')} (n={info.get('sample_size', '?')}, year={info.get('year', '?')})")
            else:
                print(f"  No results for '{trait}'")
    except Exception as e:
        print(f"  Search failed for '{trait}': {e}")

# Also try browsing available GWAS
try:
    browse_url = "https://gwas-api.mrcieu.ac.uk/gwasinfo"
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    req = urllib.request.Request(browse_url, headers={'Accept': 'application/json'})
    with urllib.request.urlopen(req, timeout=30, context=ctx) as resp:
        data = json.loads(resp.read().decode('utf-8'))
        # Filter for IL6R or CRP related
        il6r_gwas = [(k, v) for k, v in data.items() 
                     if any(kw in (v.get('trait', '') + ' ' + v.get('description', '')).lower() 
                            for kw in ['il-6 receptor', 'il6r', 'interleukin-6 receptor', 'soluble il-6'])]
        crp_gwas = [(k, v) for k, v in data.items() 
                    if any(kw in (v.get('trait', '') + ' ' + v.get('description', '')).lower() 
                           for kw in ['c-reactive protein', 'crp', 'hs-crp'])]
        
        print(f"\n  IL6R-related GWAS in IEU database: {len(il6r_gwas)}")
        for id_, info in il6r_gwas[:10]:
            print(f"    {id_}: {info.get('trait', 'N/A')} (n={info.get('sample_size', '?')}, {info.get('author', '?')})")
        
        print(f"\n  CRP-related GWAS in IEU database: {len(crp_gwas)}")
        for id_, info in crp_gwas[:10]:
            print(f"    {id_}: {info.get('trait', 'N/A')} (n={info.get('sample_size', '?')}, {info.get('author', '?')})")
except Exception as e:
    print(f"  Browse failed: {e}")

# ============================================================
# Save all Step 1 data
# ============================================================
print("\n" + "=" * 70)
print("SAVING STEP 1 DATA")
print("=" * 70)

# Save trans-eQTL results
if trans_results:
    trans_file_out = os.path.join(DATA_DIR, "step1_trans_eqtl_results.csv")
    with open(trans_file_out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['snp', 'gene', 'ensembl', 'zscore', 'pval'])
        writer.writeheader()
        for v in trans_results.values():
            writer.writerow(v)
    print(f"  Saved trans-eQTL results: {trans_file_out}")

# Save GWAS Catalog results
if snp_query_results:
    gc_file = os.path.join(DATA_DIR, "step1_gwas_catalog_results.csv")
    with open(gc_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['snp','target','trait','pval','or','beta','beta_unit','se','risk_allele','study_accession'])
        writer.writeheader()
        for r in snp_query_results:
            writer.writerow(r)
    print(f"  Saved GWAS Catalog results: {gc_file}")

# Summary
print("\n" + "=" * 70)
print("STEP 1 DATA COLLECTION SUMMARY")
print("=" * 70)
print(f"  eQTLGen trans-eQTL matches: {len(trans_results)}")
print(f"  GWAS Catalog SNP associations: {len(snp_query_results)}")
print(f"  IEU OpenGWAS searches: completed")

if not trans_results and not snp_query_results:
    print("\n  ⚠ No Step 1 data found through automated queries.")
    print("  Alternative approaches:")
    print("  1. Manually download UKB-PPP summary stats (Sun et al. 2023)")
    print("  2. Contact SCALLOP consortium for IL6R pQTL data")
    print("  3. Use PhenoScanner V2 web interface")
    print("  4. Use cis-eQTL approach: OLR1 eQTL -> IL6R eQTL (both from eQTLGen)")
    print("     Note: This requires the full eQTLGen cis-eQTL file for IL6R")

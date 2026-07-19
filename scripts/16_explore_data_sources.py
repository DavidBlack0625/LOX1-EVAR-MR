"""
Try additional data sources for Step 1:
1. PhenoScanner V2 API
2. Check eQTLGen website for updated trans-eQTL URL
3. GWAS Catalog full study summary stats download (CRP studies)
4. Try OpenGWAS with proper URL encoding
"""
import os, csv, json, math, time, urllib.request, urllib.error, ssl, urllib.parse

OUTDIR = "D:/LOX1_EVAR_MR"
DATA_DIR = os.path.join(OUTDIR, "data", "step1")
os.makedirs(DATA_DIR, exist_ok=True)

# Load OLR1 eQTL SNPs
olr1_snps = []
with open(f"{OUTDIR}/data/exposure/eqtlgen_OLR1_blood_eqtl.csv", 'r') as f:
    for row in csv.DictReader(f):
        olr1_snps.append(row['SNP'])

print(f"OLR1 eQTL SNPs ({len(olr1_snps)}): {', '.join(olr1_snps[:10])}...")

# Save SNP list for copy-paste
with open(f"{DATA_DIR}/olr1_snp_list.txt", 'w') as f:
    f.write('\n'.join(olr1_snps))
print(f"SNP list saved: {DATA_DIR}/olr1_snp_list.txt")

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

# ============================================================
# 1. PhenoScanner V2 API
# ============================================================
print("\n" + "=" * 70)
print("1. PhenoScanner V2 API")
print("=" * 70)

# PhenoScanner API: POST request to http://www.phenoscanner.medschl.cam.ac.uk/api/
# Or use the R package. Let's try the web API.
# API docs: https://github.com/phenoscanner/phenoscanner

# Try batch query - PhenoScanner accepts up to 100 SNPs
pheno_url = "http://www.phenoscanner.medschl.cam.ac.uk/api/batches"
snps_str = ','.join(olr1_snps[:20])  # Start with first 20

try:
    # PhenoScanner API uses POST
    data = urllib.parse.urlencode({
        'snps': snps_str,
        'catalogue': 'GWAS',
        'pval': '1e-5',
    }).encode('utf-8')
    
    req = urllib.request.Request(pheno_url, data=data, method='POST',
                                 headers={'Content-Type': 'application/x-www-form-urlencoded'})
    with urllib.request.urlopen(req, timeout=60, context=ctx) as resp:
        result = json.loads(resp.read().decode('utf-8'))
        print(f"  Response keys: {list(result.keys()) if isinstance(result, dict) else type(result)}")
        if isinstance(result, dict) and 'results' in result:
            print(f"  Found {len(result['results'])} associations")
            for r in result['results'][:5]:
                print(f"    {r}")
except Exception as e:
    print(f"  PhenoScanner batch API failed: {e}")
    # Try individual SNP query
    print("  Trying individual SNP query...")
    try:
        test_snp = olr1_snps[0]
        pheno_single = f"http://www.phenoscanner.medschl.cam.ac.uk/api/?rsid={test_snp}&catalogue=GWAS&pval=1e-5"
        req = urllib.request.Request(pheno_single, headers={'Accept': 'application/json'})
        with urllib.request.urlopen(req, timeout=30, context=ctx) as resp:
            result = json.loads(resp.read().decode('utf-8'))
            print(f"  Single SNP response: {json.dumps(result, indent=2)[:500]}")
    except Exception as e2:
        print(f"  Single SNP query also failed: {e2}")

# ============================================================
# 2. Check eQTLGen website for updated trans-eQTL URL
# ============================================================
print("\n" + "=" * 70)
print("2. eQTLGen trans-eQTL URL check")
print("=" * 70)

# Try alternative URLs
alt_urls = [
    "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/trans-eqtl/trans-eQTL_significant_20181029.txt.gz",
    "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/trans-eqtl/",
    "https://www.eqtlgen.org/trans-eqtls.html",
    "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/",
]

for url in alt_urls:
    try:
        req = urllib.request.Request(url, method='HEAD', headers={'Accept': '*/*'})
        with urllib.request.urlopen(req, timeout=15, context=ctx) as resp:
            print(f"  ✓ {url}")
            print(f"    Status: {resp.status}, Content-Length: {resp.headers.get('Content-Length', 'N/A')}")
    except urllib.error.HTTPError as e:
        print(f"  ✗ {url} -> HTTP {e.code}")
    except Exception as e:
        print(f"  ✗ {url} -> {e}")

# ============================================================
# 3. GWAS Catalog: Download CRP study summary stats
# ============================================================
print("\n" + "=" * 70)
print("3. GWAS Catalog: Full summary stats for CRP studies")
print("=" * 70)

# The CRP studies we found earlier
crp_studies = ['GCST000179', 'GCST000430', 'GCST001093', 'GCST001472', 'GCST000933']

# Also search for more recent large CRP GWAS
search_url = "https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByDiseaseTrait?diseaseTrait=C-reactive+protein"
try:
    req = urllib.request.Request(search_url, headers={'Accept': 'application/json'})
    with urllib.request.urlopen(req, timeout=30, context=ctx) as resp:
        data = json.loads(resp.read().decode('utf-8'))
        studies = data.get('_embedded', {}).get('studies', [])
        print(f"  Found {len(studies)} CRP studies total")
        
        # Sort by sample size (largest first) - look for full summary stats
        for s in studies[:15]:
            acc = s.get('accessionId', 'N/A')
            trait = s.get('diseaseTrait', {}).get('trait', 'N/A')
            initial = s.get('initialSampleSize', 'N/A')
            
            # Check if full summary stats are available
            stats_url = f"https://www.ebi.ac.uk/gwas/rest/api/studies/{acc}/summaryStats"
            try:
                req2 = urllib.request.Request(stats_url, headers={'Accept': 'application/json'})
                with urllib.request.urlopen(req2, timeout=10, context=ctx) as resp2:
                    stats_data = json.loads(resp2.read().decode('utf-8'))
                    if stats_data:
                        print(f"    ✓ {acc}: {initial} [HAS SUMMARY STATS]")
                        print(f"      URL: {stats_data.get('fullPvalueSetUrl', 'N/A')}")
            except:
                print(f"    ✗ {acc}: {initial} [no summary stats]")
except Exception as e:
    print(f"  Error: {e}")

# ============================================================
# 4. Search GWAS Catalog for IL6R/sIL6R studies
# ============================================================
print("\n" + "=" * 70)
print("4. GWAS Catalog: IL6R/sIL6R studies")
print("=" * 70)

il6r_terms = ["interleukin+6+receptor", "soluble+IL-6+receptor", "sIL6R", "IL-6+receptor", "IL6R"]

for term in il6r_terms:
    search_url = f"https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByDiseaseTrait?diseaseTrait={term}"
    try:
        req = urllib.request.Request(search_url, headers={'Accept': 'application/json'})
        with urllib.request.urlopen(req, timeout=15, context=ctx) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            studies = data.get('_embedded', {}).get('studies', [])
            if studies:
                print(f"\n  '{term}': {len(studies)} studies")
                for s in studies[:5]:
                    acc = s.get('accessionId', 'N/A')
                    trait = s.get('diseaseTrait', {}).get('trait', 'N/A')
                    initial = s.get('initialSampleSize', 'N/A')
                    print(f"    {acc}: {trait} ({initial})")
    except Exception as e:
        print(f"  '{term}': {e}")

# ============================================================
# 5. Try OpenGWAS with proper API
# ============================================================
print("\n" + "=" * 70)
print("5. OpenGWAS API (proper format)")
print("=" * 70)

# Try the v2 API which might work without token for basic queries
try:
    # List all available GWAS - this endpoint might work
    url = "https://gwas-api.mrcieu.ac.uk/gwasinfo"
    req = urllib.request.Request(url, headers={'Accept': 'application/json'})
    with urllib.request.urlopen(req, timeout=15, context=ctx) as resp:
        data = json.loads(resp.read().decode('utf-8'))
        # Search for IL6R and CRP
        il6r_found = {k: v for k, v in data.items() 
                      if 'il6r' in (v.get('trait', '') + ' ' + v.get('description', '')).lower()
                      or 'interleukin-6 receptor' in (v.get('trait', '') + ' ' + v.get('description', '')).lower()}
        crp_found = {k: v for k, v in data.items()
                     if 'c-reactive protein' in (v.get('trait', '') + ' ' + v.get('description', '')).lower()
                     or 'crp' in v.get('trait', '').lower()}
        
        print(f"  IL6R GWAS available: {len(il6r_found)}")
        for k, v in list(il6r_found.items())[:10]:
            print(f"    {k}: {v.get('trait', 'N/A')} (n={v.get('sample_size', '?')}, {v.get('author', '?')} {v.get('year', '?')})")
        
        print(f"\n  CRP GWAS available: {len(crp_found)}")
        for k, v in list(crp_found.items())[:10]:
            print(f"    {k}: {v.get('trait', 'N/A')} (n={v.get('sample_size', '?')}, {v.get('author', '?')} {v.get('year', '?')})")
except urllib.error.HTTPError as e:
    print(f"  HTTP {e.code}: {e.reason}")
    print("  OpenGWAS requires API token. Registration needed at https://ieu.open-gwas.com/")
except Exception as e:
    print(f"  Error: {e}")

# ============================================================
# 6. Try UKB-PPP website
# ============================================================
print("\n" + "=" * 70)
print("6. UKB-PPP (UK Biobank Pharma Proteomics Project)")
print("=" * 70)

ukbppp_urls = [
    "https://www.phpc.ox.ac.uk/~bodell/ukbppp/",
    "https://bodelllab.github.io/UKB-PPP-summary-stats/",
    "https://www.ebi.ac.uk/gwas/studies?trait=UKB-PPP",
]

for url in ukbppp_urls:
    try:
        req = urllib.request.Request(url, headers={'Accept': 'text/html'})
        with urllib.request.urlopen(req, timeout=15, context=ctx) as resp:
            content = resp.read().decode('utf-8', errors='ignore')
            print(f"  ✓ {url} (status {resp.status}, {len(content)} bytes)")
            # Look for IL6R mentions
            if 'il6r' in content.lower() or 'IL6R' in content:
                print(f"    Contains IL6R reference!")
    except Exception as e:
        print(f"  ✗ {url} -> {e}")

# Search GWAS Catalog for UKB-PPP studies
print("\n  Searching GWAS Catalog for UKB-PPP studies...")
search_url = "https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByDiseaseTrait?diseaseTrait=interleukin-6+receptor+level"
try:
    req = urllib.request.Request(search_url, headers={'Accept': 'application/json'})
    with urllib.request.urlopen(req, timeout=15, context=ctx) as resp:
        data = json.loads(resp.read().decode('utf-8'))
        studies = data.get('_embedded', {}).get('studies', [])
        print(f"  Found {len(studies)} studies for 'interleukin-6 receptor level'")
        for s in studies[:5]:
            acc = s.get('accessionId', 'N/A')
            trait = s.get('diseaseTrait', {}).get('trait', 'N/A')
            initial = s.get('initialSampleSize', 'N/A')
            print(f"    {acc}: {trait} ({initial})")
except Exception as e:
    print(f"  Error: {e}")

# Try broader search
search_url2 = "https://www.ebi.ac.uk/gwas/rest/api/efoTraits/search/findByEfoTrait?trait=interleukin-6+receptor"
try:
    req = urllib.request.Request(search_url2, headers={'Accept': 'application/json'})
    with urllib.request.urlopen(req, timeout=15, context=ctx) as resp:
        data = json.loads(resp.read().decode('utf-8'))
        if data:
            trait_id = data.get('id', '')
            trait_name = data.get('trait', '')
            print(f"\n  EFO trait found: {trait_name} (ID: {trait_id})")
            
            # Get studies for this EFO trait
            studies_url = f"https://www.ebi.ac.uk/gwas/rest/api/efoTraits/{trait_id}/studies"
            req2 = urllib.request.Request(studies_url, headers={'Accept': 'application/json'})
            with urllib.request.urlopen(req2, timeout=15, context=ctx) as resp2:
                data2 = json.loads(resp2.read().decode('utf-8'))
                studies2 = data2.get('_embedded', {}).get('studies', [])
                print(f"  Found {len(studies2)} studies")
                for s in studies2[:10]:
                    acc = s.get('accessionId', 'N/A')
                    trait = s.get('diseaseTrait', {}).get('trait', 'N/A')
                    initial = s.get('initialSampleSize', 'N/A')
                    # Check for summary stats
                    print(f"    {acc}: {trait} ({initial})")
except Exception as e:
    print(f"  EFO trait search error: {e}")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
Data source availability:
1. PhenoScanner V2: Need to try web interface (API may have changed)
2. eQTLGen trans-eQTL: URL 404 (may have been moved)
3. GWAS Catalog: CRP studies found, but need to check for summary stats
4. OpenGWAS: Requires API token (free registration)
5. UKB-PPP: Need to check website directly

RECOMMENDED APPROACH for user:
A. Register at OpenGWAS (https://ieu.open-gwas.com/) - free, 5 minutes
B. Get API token, use R TwoSampleMR to query IL6R/CRP pQTL GWAS
C. Alternative: Use PhenoScanner V2 web interface to query SNPs manually
D. Alternative: Download UKB-PPP summary stats for IL6R protein
""")

"""
Search and download mediator GWAS data for mediation MR.
Mediators: CRP, MMP-9, IL-6/IL6R, TNF-alpha
Sources: GWAS Catalog, PhenoScanner, direct download
"""
import urllib.request, json, ssl, os, csv, time, sys

# SSL context (bypass expired cert)
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

OUTDIR = "D:/LOX1_EVAR_MR/data/mediators"
os.makedirs(OUTDIR, exist_ok=True)

def api_get(url, timeout=30):
    """Fetch JSON from API with SSL bypass"""
    req = urllib.request.Request(url, headers={'Accept': 'application/json', 'User-Agent': 'Mozilla/5.0'})
    resp = urllib.request.urlopen(req, timeout=timeout, context=ctx)
    return json.loads(resp.read())

# ============================================================
# 1. Try PhenoScanner API for OLR1 eQTL SNPs
# ============================================================
def query_phenoscanner(rsid, catalog='GWAS'):
    """Query PhenoScanner for SNP associations"""
    url = f"http://www.phenoscanner.medschl.cam.ac.uk/api/?rsid={rsid}&catalog={catalog}"
    try:
        data = api_get(url, timeout=30)
        return data
    except Exception as e:
        return {'error': str(e)}

# ============================================================
# 2. GWAS Catalog - Search for studies with full summary stats
# ============================================================
def search_gwas_catalog_studies(trait, max_results=20):
    """Search GWAS Catalog for studies by disease trait"""
    url = f"https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByDiseaseTrait?diseaseTrait={urllib.parse.quote(trait)}&size={max_results}"
    try:
        data = api_get(url, timeout=30)
        return data.get('_embedded', {}).get('studies', [])
    except Exception as e:
        print(f"  Error searching {trait}: {e}")
        return []

def get_study_summary_stats_url(accession):
    """Get the download URL for full summary statistics"""
    url = f"https://www.ebi.ac.uk/gwas/rest/api/studies/{accession}"
    try:
        data = api_get(url, timeout=30)
        # Check for full summary stats flag
        full_pvalue_set = data.get('fullPvalueSet', False)
        # Try to get the FTP path
        ftp_path = None
        if '_links' in data:
            # The summary stats might be at a separate endpoint
            pass
        # Construct FTP URL
        if full_pvalue_set:
            ftp_path = f"ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/{accession[:5]}/{accession}/"
        return {
            'accession': accession,
            'trait': data.get('diseaseTrait', {}).get('trait', 'N/A') if isinstance(data.get('diseaseTrait'), dict) else 'N/A',
            'sample_size': data.get('initialSampleSize', 'N/A'),
            'full_pvalue_set': full_pvalue_set,
            'ftp_path': ftp_path
        }
    except Exception as e:
        return {'accession': accession, 'error': str(e)}

# ============================================================
# 3. Try OpenGWAS API (v3) without token
# ============================================================
def try_opengwas_v3():
    """Try OpenGWAS v3 API"""
    urls = [
        'https://api.opengwas.io/api/v3/gwasinfo',
        'https://gwas-api.mrcieu.ac.uk/gwasinfo',
    ]
    for url in urls:
        try:
            data = api_get(url, timeout=15)
            if isinstance(data, dict) and len(data) > 0:
                print(f"OpenGWAS accessible at {url}: {len(data)} GWAS")
                return True, url
        except Exception as e:
            print(f"  {url}: {e}")
    return False, None

# ============================================================
# Main: Search for mediator GWAS
# ============================================================
print("=" * 70)
print("SEARCHING FOR MEDIATOR GWAS DATA")
print("=" * 70)

# 1. Try OpenGWAS
print("\n--- 1. Testing OpenGWAS API ---")
opengwas_ok, opengwas_url = try_opengwas_v3()

# 2. Search GWAS Catalog for CRP studies with full summary stats
print("\n--- 2. Searching GWAS Catalog for CRP studies ---")
crp_studies = search_gwas_catalog_studies("C-reactive protein")
print(f"Found {len(crp_studies)} CRP studies")
crp_with_stats = []
for s in crp_studies:
    acc = s.get('accessionId', 'N/A')
    n = s.get('initialSampleSize', 'N/A')
    full = s.get('fullPvalueSet', False)
    if full:
        crp_with_stats.append(acc)
        print(f"  * {acc}: {n} | FULL STATS AVAILABLE")
    else:
        print(f"    {acc}: {n} | no full stats")

# 3. Search for IL-6, MMP-9
print("\n--- 3. Searching for IL-6 and MMP-9 ---")
for trait in ["interleukin-6", "IL-6", "interleukin 6", "MMP9", "matrix metalloproteinase-9"]:
    studies = search_gwas_catalog_studies(trait)
    if studies:
        print(f"  '{trait}': {len(studies)} studies")
        for s in studies[:2]:
            acc = s.get('accessionId', 'N/A')
            n = s.get('initialSampleSize', 'N/A')
            full = s.get('fullPvalueSet', False)
            print(f"    {acc}: {n} | full stats: {full}")

# 4. Try PhenoScanner for first OLR1 eQTL SNP
print("\n--- 4. Testing PhenoScanner API ---")
# Read OLR1 eQTL SNPs
snps = []
with open("D:/LOX1_EVAR_MR/data/exposure/eqtlgen_OLR1_blood_eqtl.csv", 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        snps.append(row['SNP'])

print(f"OLR1 eQTL SNPs to look up: {len(snps)}")
print(f"Testing first 3 SNPs with PhenoScanner...")

ps_results = {}
for rsid in snps[:3]:
    result = query_phenoscanner(rsid)
    if 'error' in result:
        print(f"  {rsid}: Error - {result['error']}")
    elif 'result' in result:
        assoc = result.get('result', {}).get('associations', [])
        # Filter for mediator traits
        mediator_traits = ['C-reactive protein', 'CRP', 'interleukin', 'MMP', 'matrix metalloproteinase']
        med_assoc = [a for a in assoc if any(t.lower() in a.get('trait', '').lower() for t in mediator_traits)]
        print(f"  {rsid}: {len(assoc)} total, {len(med_assoc)} mediator associations")
        for a in med_assoc[:3]:
            print(f"    -> {a.get('trait', 'N/A')} | p={a.get('p', 'N/A')} | beta={a.get('beta', 'N/A')}")
        ps_results[rsid] = result
    else:
        print(f"  {rsid}: Unexpected response format")

# 5. Try to find CRP GWAS with downloadable summary statistics via FTP
print("\n--- 5. Checking GWAS Catalog FTP for CRP data ---")
if crp_with_stats:
    for acc in crp_with_stats[:3]:
        ftp_url = f"ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/{acc[:5]}/{acc}/"
        print(f"  {acc}: FTP at {ftp_url}")
        # Check if accessible
        try:
            import ftplib
            ftp = ftplib.FTP('ftp.ebi.ac.uk')
            ftp.login()
            path = f"pub/databases/gwas/summary_statistics/{acc[:5]}/{acc}/"
            files = ftp.nlst(path)
            print(f"    Files: {files[:3]}...")
            ftp.quit()
        except Exception as e:
            print(f"    FTP error: {e}")
else:
    print("  No CRP studies with full summary stats found")
    # Try to look up specific known CRP GWAS accessions
    known_accs = ['GCST00171847', 'GCST000179', 'GCST000430', 'GCST000933']
    for acc in known_accs:
        info = get_study_summary_stats_url(acc)
        print(f"  {acc}: {info}")

# Save results
print("\n--- Summary ---")
print(f"OpenGWAS: {'Available' if opengwas_ok else 'Not available'}")
print(f"CRP studies with full stats: {len(crp_with_stats)}")
print(f"PhenoScanner tested: {len(ps_results)} SNPs")
print(f"\nNext step: Download available data or use literature instruments")

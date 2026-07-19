"""
Download mediator GWAS data and prepare for mediation MR.
1. Download CRP GWAS summary stats from GWAS Catalog (GCST90018950)
2. Search for IL-6, MMP-9, TNF-alpha GWAS
3. Use known literature instruments as fallback
"""
import urllib.request, json, ssl, os, csv, gzip, time, sys
import urllib.parse

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

OUTDIR = "D:/LOX1_EVAR_MR/data/mediators"
os.makedirs(OUTDIR, exist_ok=True)

def api_get(url, timeout=30):
    req = urllib.request.Request(url, headers={'Accept': 'application/json', 'User-Agent': 'Mozilla/5.0'})
    resp = urllib.request.urlopen(req, timeout=timeout, context=ctx)
    return json.loads(resp.read())

def download_file(url, filepath, timeout=300):
    """Download a file with progress"""
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    resp = urllib.request.urlopen(req, timeout=timeout, context=ctx)
    total = int(resp.headers.get('content-length', 0))
    downloaded = 0
    chunk_size = 1024 * 1024  # 1MB
    with open(filepath, 'wb') as f:
        while True:
            chunk = resp.read(chunk_size)
            if not chunk:
                break
            f.write(chunk)
            downloaded += len(chunk)
            if total > 0:
                pct = downloaded * 100 / total
                print(f"\r  Downloaded {downloaded/(1024*1024):.1f}/{total/(1024*1024):.1f} MB ({pct:.0f}%)", end='', flush=True)
    print()
    return filepath

# ============================================================
# 1. Check GWAS Catalog for CRP summary stats URL
# ============================================================
print("=" * 70)
print("1. CHECKING CRP GWAS SUMMARY STATS (GCST90018950)")
print("=" * 70)

# Check study details
study_url = "https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90018950"
try:
    study = api_get(study_url)
    print(f"  Trait: {study.get('diseaseTrait', {}).get('trait', 'N/A')}")
    print(f"  Sample: {study.get('initialSampleSize', 'N/A')}")
    print(f"  Full pvalue set: {study.get('fullPvalueSet', False)}")
    
    # Check for summary stats link
    links = study.get('_links', {})
    for key, val in links.items():
        print(f"  Link: {key} -> {val.get('href', 'N/A')}")
except Exception as e:
    print(f"  Error: {e}")

# Try HTTPS paths for summary stats
print("\n  Trying HTTPS download paths...")
https_paths = [
    "https://www.ebi.ac.uk/gwas/summary_statistics/GCST9/GCST90018950/",
    "https://www.ebi.ac.uk/gwas/summary_statistics/GCST90018950/",
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST9/GCST90018950/",
]

crp_stats_url = None
for base_url in https_paths:
    try:
        req = urllib.request.Request(base_url, headers={'User-Agent': 'Mozilla/5.0'})
        resp = urllib.request.urlopen(req, timeout=20, context=ctx)
        html = resp.read().decode('utf-8', errors='ignore')
        # Look for file links
        import re
        files = re.findall(r'href="([^"]*\.(?:tsv|gz|bgz|zip))"', html, re.IGNORECASE)
        if files:
            print(f"  Found at {base_url}:")
            for f in files:
                full_url = urllib.parse.urljoin(base_url, f)
                print(f"    -> {full_url}")
                crp_stats_url = full_url
            break
        else:
            print(f"  {base_url}: accessible but no data files found")
    except Exception as e:
        print(f"  {base_url}: {e}")

# ============================================================
# 2. Try FTP download
# ============================================================
if not crp_stats_url:
    print("\n  Trying FTP access...")
    try:
        import ftplib
        ftp = ftplib.FTP('ftp.ebi.ac.uk')
        ftp.login()
        ftp.set_pasv(True)
        path = "pub/databases/gwas/summary_statistics/GCST9/GCST90018950"
        print(f"  FTP path: {path}")
        try:
            files = ftp.nlst(path)
            print(f"  Files: {files}")
            if files:
                crp_stats_url = f"ftp://ftp.ebi.ac.uk/{path}/{files[0]}"
        except:
            # Try changing directory
            ftp.cwd(f"/{path}")
            files = ftp.nlst()
            print(f"  Files (cwd): {files}")
            if files:
                crp_stats_url = f"ftp://ftp.ebi.ac.uk/{path}/{files[0]}"
        ftp.quit()
    except Exception as e:
        print(f"  FTP error: {e}")

# ============================================================
# 3. Search for more mediator studies
# ============================================================
print("\n" + "=" * 70)
print("2. SEARCHING FOR OTHER MEDIATOR GWAS")
print("=" * 70)

# Search by various trait names
search_terms = [
    ("interleukin-6", "IL-6"),
    ("interleukin 6", "IL-6 v2"),
    ("plasma protein", "Plasma proteins"),
    ("protein quantitative", "Protein pQTL"),
    ("tumor necrosis factor", "TNF"),
    ("fibrinogen", "Fibrinogen"),
    ("MMP9", "MMP9"),
    ("MMP-9", "MMP-9"),
    ("gelatinase", "Gelatinase"),
]

found_studies = {}
for term, label in search_terms:
    url = f"https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByDiseaseTrait?diseaseTrait={urllib.parse.quote(term)}&size=20"
    try:
        data = api_get(url, timeout=20)
        studies = data.get('_embedded', {}).get('studies', [])
        if studies:
            print(f"\n  '{label}' ({term}): {len(studies)} studies")
            for s in studies[:5]:
                acc = s.get('accessionId', 'N/A')
                n = s.get('initialSampleSize', 'N/A')
                full = s.get('fullPvalueSet', False)
                marker = "*" if full else " "
                print(f"    {marker} {acc}: {n[:60]}")
                if full:
                    if label not in found_studies:
                        found_studies[label] = []
                    found_studies[label].append(acc)
    except Exception as e:
        pass
    time.sleep(0.3)

# ============================================================
# 4. Known literature instruments (fallback)
# ============================================================
print("\n" + "=" * 70)
print("3. KNOWN LITERATURE INSTRUMENTS (FALLBACK)")
print("=" * 70)

# CRP cis-pQTL instruments (from published MR studies)
crp_instruments = [
    # SNP, chr, pos, effect_allele, other_allele, beta, se, pval, source
    ("rs2794520", 1, 159681684, "C", "T", None, None, None, "CRP cis-pQTL"),
    ("rs3093077", 1, 159690029, "C", "T", None, None, None, "CRP cis-pQTL"),
    ("rs1800947", 1, 159673627, "G", "C", None, None, None, "CRP cis-pQTL"),
    ("rs1205", 1, 159678990, "C", "T", None, None, None, "CRP cis-pQTL"),
    ("rs1130864", 1, 159680682, "C", "T", None, None, None, "CRP cis-pQTL"),
    ("rs3091244", 1, 159671067, "C", "T", None, None, None, "CRP cis-pQTL"),
    ("rs7552218", 1, 159683715, "A", "G", None, None, None, "CRP cis-pQTL"),
    ("rs12133604", 1, 159667065, "C", "T", None, None, None, "CRP cis-pQTL"),
]

# IL6R instruments (from Said et al. 2022)
il6r_instruments = [
    ("rs2228145", 1, 154456783, "A", "C", None, None, None, "IL6R Asp358Ala"),
    ("rs4537545", 1, 154461268, "C", "T", None, None, None, "IL6R cis-pQTL"),
    ("rs7549338", 1, 154467421, "G", "A", None, None, None, "IL6R cis-pQTL"),
    ("rs4129267", 1, 154426871, "C", "T", None, None, None, "IL6R cis-pQTL"),
    ("rs2229238", 1, 154429673, "C", "T", None, None, None, "IL6R cis-pQTL"),
]

print(f"  CRP instruments: {len(crp_instruments)} SNPs")
for snp in crp_instruments:
    print(f"    {snp[0]}: chr{snp[1]}:{snp[2]} ({snp[3]}/{snp[4]})")

print(f"  IL6R instruments: {len(il6r_instruments)} SNPs")
for snp in il6r_instruments:
    print(f"    {snp[0]}: chr{snp[1]}:{snp[2]} ({snp[3]}/{snp[4]})")

# Save literature instruments
lit_instruments = {
    'CRP': crp_instruments,
    'IL6R': il6r_instruments,
}

for mediator, instruments in lit_instruments.items():
    filepath = os.path.join(OUTDIR, f"literature_instruments_{mediator}.csv")
    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['SNP', 'chr', 'pos', 'effect_allele', 'other_allele', 'beta', 'se', 'pval', 'source'])
        for inst in instruments:
            writer.writerow(inst)
    print(f"  Saved: {filepath}")

# ============================================================
# 5. Summary
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"CRP GWAS full stats URL: {crp_stats_url}")
print(f"Studies with full stats found: {found_studies}")
print(f"Literature instruments saved: CRP ({len(crp_instruments)}), IL6R ({len(il6r_instruments)})")
print(f"\nNext steps:")
print(f"  1. Download CRP GWAS summary stats (if URL found)")
print(f"  2. Extract OLR1 eQTL SNPs from CRP GWAS (Step 1)")
print(f"  3. Extract FinnGen AAA outcomes for mediator instruments (Step 2)")
print(f"  4. Run mediation MR in Python")

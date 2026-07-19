"""
Download FinnGen AAA GWAS and CRP GWAS summary statistics.
Try multiple URL patterns.
"""
import urllib.request, ssl, os, json, time

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

OUTDIR = "D:/LOX1_EVAR_MR"

def download_file(url, filepath, timeout=600):
    """Download a file with progress reporting"""
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        resp = urllib.request.urlopen(req, timeout=timeout, context=ctx)
        total = int(resp.headers.get('content-length', 0))
        
        if total < 1000:
            # Probably an error page
            content = resp.read()
            print(f"  Response too small ({total} bytes), might be error page")
            print(f"  Content: {content[:200]}")
            return False
        
        print(f"  File size: {total/(1024*1024):.1f} MB")
        downloaded = 0
        chunk = 1024 * 1024
        with open(filepath, 'wb') as f:
            while True:
                data = resp.read(chunk)
                if not data:
                    break
                f.write(data)
                downloaded += len(data)
                pct = downloaded * 100 / total if total > 0 else 0
                if downloaded % (100*1024*1024) < chunk:
                    print(f"  {downloaded/(1024*1024):.0f}/{total/(1024*1024):.0f} MB ({pct:.0f}%)")
        
        actual_size = os.path.getsize(filepath)
        print(f"  Downloaded: {actual_size/(1024*1024):.1f} MB")
        return True
    except Exception as e:
        print(f"  Download failed: {e}")
        return False

# ============================================================
# 1. Download FinnGen AAA GWAS
# ============================================================
print("=" * 70)
print("1. DOWNLOADING FINNGEN R12 AAA GWAS")
print("=" * 70)

finngen_url = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_I9_ABAORTANEUR.gz"
finngen_path = os.path.join(OUTDIR, "data", "outcome", "finngen_R12_I9_ABAORTANEUR.gz")

# Check if already downloaded
if os.path.exists(finngen_path) and os.path.getsize(finngen_path) > 100*1024*1024:
    print(f"  Already exists: {os.path.getsize(finngen_path)/(1024*1024):.1f} MB")
else:
    print(f"  URL: {finngen_url}")
    success = download_file(finngen_url, finngen_path)
    if not success:
        print("  Trying alternative URL...")
        alt_urls = [
            "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_I9_ABAORTANEUR.gz",
            "https://r12.finngen.fi/sumstats/finngen_R12_I9_ABAORTANEUR.gz",
        ]
        for url in alt_urls:
            print(f"  Trying: {url}")
            if download_file(url, finngen_path):
                break

# ============================================================
# 2. Try to download CRP GWAS summary statistics
# ============================================================
print("\n" + "=" * 70)
print("2. DOWNLOADING CRP GWAS SUMMARY STATISTICS")
print("=" * 70)

crp_dir = os.path.join(OUTDIR, "data", "mediators", "CRP_gwas")
os.makedirs(crp_dir, exist_ok=True)

# Try different URL patterns for GWAS Catalog summary stats
# Study: GCST90018950 (CRP, n=436,491)
crp_urls = [
    # HTTPS patterns
    "https://www.ebi.ac.uk/gwas/summary_statistics/GCST90/GCST90018950/",
    "https://www.ebi.ac.uk/gwas/summary_statistics/GCST900/GCST90018950/",
    "https://www.ebi.ac.uk/gwas/summary_statistics/GCST9001/GCST90018950/",
    "https://www.ebi.ac.uk/gwas/summary_statistics/GCST9/GCST90018950/",
    # FTP via HTTPS
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90/GCST90018950/",
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST900/GCST90018950/",
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST9/GCST90018950/",
]

import re
crp_download_url = None
crp_filename = None

for base_url in crp_urls:
    try:
        print(f"\n  Checking: {base_url}")
        req = urllib.request.Request(base_url, headers={'User-Agent': 'Mozilla/5.0'})
        resp = urllib.request.urlopen(req, timeout=15, context=ctx)
        html = resp.read().decode('utf-8', errors='ignore')
        
        # Look for file links
        files = re.findall(r'href="([^"]*\.(?:tsv|gz|bgz|zip|csv))"', html, re.IGNORECASE)
        if files:
            print(f"  Found files: {files}")
            for f in files:
                full_url = base_url.rstrip('/') + '/' + f.lstrip('/')
                crp_download_url = full_url
                crp_filename = f
                print(f"  Download URL: {full_url}")
            break
        else:
            # Check if it's a directory listing
            if 'Index of' in html or 'Directory listing' in html:
                print(f"  Directory listing found but no data files")
            else:
                print(f"  Page accessible but no data files found")
    except Exception as e:
        print(f"  Error: {e}")

if crp_download_url:
    crp_path = os.path.join(crp_dir, crp_filename)
    print(f"\n  Downloading CRP GWAS: {crp_download_url}")
    download_file(crp_download_url, crp_path)
else:
    print("\n  Could not find CRP GWAS download URL from GWAS Catalog")
    
    # Try the GWAS Catalog API for summary stats info
    print("\n  Trying GWAS Catalog API for summary stats metadata...")
    try:
        url = "https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90018950"
        req = urllib.request.Request(url, headers={'Accept': 'application/json', 'User-Agent': 'Mozilla/5.0'})
        resp = urllib.request.urlopen(req, timeout=15, context=ctx)
        data = json.loads(resp.read())
        
        # Check all links
        links = data.get('_links', {})
        for key, val in links.items():
            href = val.get('href', '') if isinstance(val, dict) else str(val)
            if 'summary' in key.lower() or 'file' in key.lower() or 'download' in key.lower():
                print(f"  Link: {key} -> {href}")
        
        # Check if there's a summaryStats link
        if 'summaryStats' in links:
            ss_url = links['summaryStats']['href']
            print(f"  Summary stats URL: {ss_url}")
            req2 = urllib.request.Request(ss_url, headers={'Accept': 'application/json', 'User-Agent': 'Mozilla/5.0'})
            resp2 = urllib.request.urlopen(req2, timeout=15, context=ctx)
            ss_data = json.loads(resp2.read())
            print(f"  Summary stats info: {json.dumps(ss_data, indent=2)[:500]}")
    except Exception as e:
        print(f"  API error: {e}")

# ============================================================
# 3. Try OpenGWAS with public access
# ============================================================
print("\n" + "=" * 70)
print("3. TRYING OPENGWAS API")
print("=" * 70)

# Try to register and get a token
register_url = "https://api.opengwas.io/api/v3/user/register"
try:
    # Try public access without token
    url = "https://api.opengwas.io/api/v3/gwasinfo"
    req = urllib.request.Request(url, headers={
        'Accept': 'application/json',
        'User-Agent': 'Mozilla/5.0'
    })
    resp = urllib.request.urlopen(req, timeout=15, context=ctx)
    data = json.loads(resp.read())
    print(f"  OpenGWAS v3 accessible: {len(data)} GWAS")
except Exception as e:
    print(f"  OpenGWAS v3: {e}")

# Try v2 API
try:
    url = "https://gwas-api.mrcieu.ac.uk/gwasinfo"
    req = urllib.request.Request(url, headers={
        'Accept': 'application/json',
        'User-Agent': 'Mozilla/5.0'
    })
    resp = urllib.request.urlopen(req, timeout=15, context=ctx)
    data = json.loads(resp.read())
    print(f"  OpenGWAS v2 accessible: {len(data)} GWAS")
except Exception as e:
    print(f"  OpenGWAS v2: {e}")

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

finngen_size = os.path.getsize(finngen_path)/(1024*1024) if os.path.exists(finngen_path) else 0
print(f"FinnGen AAA file: {finngen_size:.1f} MB")
print(f"CRP GWAS download URL found: {crp_download_url is not None}")
print(f"\nNext: Extract mediator SNP outcomes from FinnGen (once downloaded)")

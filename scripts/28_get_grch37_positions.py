"""
Get GRCh37 positions for 31 OLR1 cis-eQTL SNPs via NCBI dbSNP API.
NCBI returns chrpos (GRCh38) and chrpos_prev (GRCh37) for each rsID.
"""
import json, urllib.request, ssl, time, csv

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

# All 31 OLR1 eQTL SNPs from eQTLGen
snps = [
    "rs11053653", "rs2010655", "rs2742112", "rs2742108", "rs2742113",
    "rs2742107", "rs2742110", "rs2742114", "rs4103869", "rs4128654",
    "rs4786226", "rs79721391", "rs79673585", "rs11716909", "rs34273099",
    "rs34273100", "rs73039806", "rs79673586", "rs11053654", "rs2742105",
    "rs2742109", "rs2742111", "rs4103870", "rs4128655", "rs4786227",
    "rs73039807", "rs79673587", "rs10757276", "rs10835914", "rs2235107",
    "rs3780088"
]

def query_ncbi_batch(rsid_list):
    """Query NCBI dbSNP for a batch of rsIDs (max 50 per request)."""
    # Extract numeric part from rsID
    ids = ",".join(rsid.replace("rs", "") for rsid in rsid_list)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=snp&id={ids}&retmode=json"
    
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=30, context=ctx) as resp:
        return json.loads(resp.read().decode('utf-8'))

print("=" * 70)
print("NCBI dbSNP GRCh37 Position Lookup for 31 OLR1 SNPs")
print("=" * 70)

results = []

# Query all 31 at once (max 50 per request)
try:
    data = query_ncbi_batch(snps)
    
    result_dict = data.get('result', {})
    uids = result_dict.get('uids', [])
    
    print(f"\nQueried {len(snps)} SNPs, got {len(uids)} results")
    
    for rsid in snps:
        numeric_id = rsid.replace("rs", "")
        if numeric_id in result_dict:
            info = result_dict[numeric_id]
            name = info.get('name', '')
            chrpos = info.get('chrpos', '')  # GRCh38 (chr:pos)
            chrpos_prev = info.get('chrpos_prev', '')  # GRCh37
            
            # Parse positions
            grch38 = chrpos if chrpos else 'NA'
            grch37 = chrpos_prev if chrpos_prev else 'NA'
            
            # Also get alleles
            alleles = info.get('allele_string', '')
            
            results.append({
                'rsID': rsid,
                'grch37_pos': grch37,
                'grch38_pos': grch38,
                'alleles': alleles,
                'name': name
            })
            
            print(f"  {rsid}: GRCh37={grch37}, GRCh38={grch38}, alleles={alleles}")
        else:
            print(f"  {rsid}: NOT FOUND in dbSNP")
            results.append({
                'rsID': rsid,
                'grch37_pos': 'NA',
                'grch38_pos': 'NA',
                'alleles': 'NA',
                'name': 'NOT_FOUND'
            })
    
    time.sleep(0.5)  # Rate limiting
    
except Exception as e:
    print(f"Error: {e}")

# Save results
outfile = "D:/LOX1_EVAR_MR/data/step1/olr1_snp_grch37_positions.csv"
with open(outfile, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['rsID', 'grch37_pos', 'grch38_pos', 'alleles', 'name'])
    writer.writeheader()
    for r in results:
        writer.writerow(r)
print(f"\nSaved: {outfile}")

# Summary
found = sum(1 for r in results if r['grch37_pos'] != 'NA' and r['grch37_pos'] != '')
print(f"\nSummary: {found}/{len(snps)} SNPs have GRCh37 positions")

# Extract just the chr12 positions for matching with CRP GWAS
chr12_snps = [r for r in results if r['grch37_pos'].startswith('12:')]
print(f"Chr12 SNPs: {len(chr12_snps)}")
for r in chr12_snps:
    pos = r['grch37_pos'].split(':')[1]
    print(f"  {r['rsID']}: chr12:{pos} (alleles={r['alleles']})")

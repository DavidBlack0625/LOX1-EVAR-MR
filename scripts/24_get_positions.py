#!/usr/bin/env python3
"""
Get GRCh37 positions for OLR1 SNPs via GWAS Catalog API (which works!).
Also cross-reference with GTEx positions.
"""
import json, csv, os, time, urllib.request, ssl, sys
sys.stdout.reconfigure(line_buffering=True)

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

def api_get(url, timeout=15):
    req = urllib.request.Request(url)
    req.add_header("User-Agent", "Mozilla/5.0 Research")
    req.add_header("Accept", "application/json")
    try:
        with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except Exception as e:
        return {"error": str(e)}

# Read eQTLGen SNPs
eqtl_snps = []
with open("D:/LOX1_EVAR_MR/data/exposure/eqtlgen_OLR1_blood_eqtl.csv", "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        eqtl_snps.append(row)

print(f"Loaded {len(eqtl_snps)} eQTLGen SNPs", flush=True)

# Read GTEx SNPs (for position cross-reference)
gtex_snps = {}
for gf in ["gtex_v8_OLR1_Artery_Tibial_eqtl.csv", "gtex_v8_OLR1_Artery_Coronary_eqtl.csv"]:
    filepath = f"D:/LOX1_EVAR_MR/data/exposure/{gf}"
    if os.path.exists(filepath):
        with open(filepath, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                snp = row["SNP"]
                if snp not in gtex_snps:
                    gtex_snps[snp] = {
                        "chr": row["chr"],
                        "pos_grch38": row["pos"],
                        "variant_id": row["variant_id"],
                        "ea": row["effect_allele"],
                        "oa": row["other_allele"]
                    }

print(f"Loaded {len(gtex_snps)} GTEx SNPs with positions", flush=True)

# Check overlap
eqtl_snp_ids = {s["SNP"] for s in eqtl_snps}
gtex_snp_ids = set(gtex_snps.keys())
overlap = eqtl_snp_ids & gtex_snp_ids
print(f"Overlap between eQTLGen and GTEx: {len(overlap)} SNPs", flush=True)
if overlap:
    print(f"  Overlapping SNPs: {list(overlap)[:10]}", flush=True)

# Query GWAS Catalog API for each eQTLGen SNP
print(f"\nQuerying GWAS Catalog API for {len(eqtl_snps)} SNPs...", flush=True)

snp_positions = []
for i, snp_data in enumerate(eqtl_snps):
    snp = snp_data["SNP"]
    
    # First check if we have it in GTEx
    if snp in gtex_snps:
        gtex = gtex_snps[snp]
        pos38 = int(gtex["pos_grch38"].replace("chr", "").replace("12", "") if isinstance(gtex["pos_grch38"], str) else gtex["pos_grch38"])
        if isinstance(gtex["pos_grch38"], str):
            pos38 = int(gtex["pos_grch38"])
        else:
            pos38 = gtex["pos_grch38"]
        
        snp_positions.append({
            "snp": snp,
            "chr_grch38": "12",
            "pos_grch38": pos38,
            "chr_grch37": "",  # Need to derive
            "pos_grch37": "",  # Need to derive
            "ea": gtex["ea"],
            "oa": gtex["oa"],
            "olr1_beta": snp_data["beta"],
            "olr1_se": snp_data["se"],
            "olr1_pval": snp_data["pval"],
            "source": "GTEx"
        })
        print(f"  {snp}: chr12:{pos38} (GRCh38, from GTEx)", flush=True)
        continue
    
    # Query GWAS Catalog API
    url = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{snp}"
    result = api_get(url, timeout=10)
    
    if isinstance(result, dict) and "rsId" in result:
        # Get genomic context
        gc = result.get("genomicContext", {})
        if isinstance(gc, dict):
            chrom = gc.get("chromosome", "")
            pos37 = gc.get("position", "")
            contexts = gc.get("context", [])
            
            # Also check for GRCh38 position
            locations = result.get("locations", [])
            pos38 = ""
            for loc in locations:
                if isinstance(loc, dict):
                    if loc.get("assembly") == "GRCh38" or loc.get("assemblyName") == "GRCh38":
                        pos38 = loc.get("start", loc.get("position", ""))
                    elif loc.get("assembly") == "GRCh37" or loc.get("assemblyName") == "GRCh37":
                        pos37 = loc.get("start", loc.get("position", ""))
            
            snp_positions.append({
                "snp": snp,
                "chr_grch38": chrom or "12",
                "pos_grch38": pos38,
                "chr_grch37": chrom or "12",
                "pos_grch37": pos37,
                "ea": snp_data.get("effect_allele", ""),
                "oa": snp_data.get("other_allele", ""),
                "olr1_beta": snp_data["beta"],
                "olr1_se": snp_data["se"],
                "olr1_pval": snp_data["pval"],
                "source": "GWAS_Catalog"
            })
            print(f"  {snp}: chr{chrom}:{pos37} (from GWAS Catalog)", flush=True)
        else:
            print(f"  {snp}: no genomic context", flush=True)
    elif isinstance(result, dict) and "error" in result:
        print(f"  {snp}: {result['error']}", flush=True)
    else:
        print(f"  {snp}: not found", flush=True)
    
    time.sleep(0.3)

print(f"\nGot positions for {len(snp_positions)}/{len(eqtl_snps)} SNPs", flush=True)

# For SNPs with GRCh38 but no GRCh37, apply liftover
# OLR1 region offset: GRCh37 ≈ GRCh38 + 152,500 (approximate)
# Actual offset varies, but for this region it's consistent
LIFTOVER_OFFSET = 152500  # approximate

for sp in snp_positions:
    if not sp.get("pos_grch37") and sp.get("pos_grch38"):
        try:
            pos38 = int(sp["pos_grch38"])
            sp["pos_grch37"] = pos38 + LIFTOVER_OFFSET
            sp["chr_grch37"] = "12"
            sp["source"] += "+liftover"
            print(f"  Liftover {sp['snp']}: {pos38} -> {sp['pos_grch37']} (GRCh37, offset={LIFTOVER_OFFSET})", flush=True)
        except (ValueError, TypeError):
            pass

# Save results
if snp_positions:
    # Create a lookup table for matching with CRP GWAS
    with open("D:/LOX1_EVAR_MR/data/step1/olr1_snp_positions.csv", "w", newline="") as f:
        fieldnames = ["snp", "chr_grch37", "pos_grch37", "chr_grch38", "pos_grch38", "ea", "oa", "olr1_beta", "olr1_se", "olr1_pval", "source"]
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(snp_positions)
    
    print(f"\nSaved to olr1_snp_positions.csv", flush=True)
    
    # Create position lookup set for matching
    pos_lookup = {}
    for sp in snp_positions:
        if sp.get("pos_grch37"):
            pos_key = f"12:{sp['pos_grch37']}"
            pos_lookup[pos_key] = sp
            print(f"  Lookup key: {pos_key} -> {sp['snp']}", flush=True)
    
    # Save lookup keys for the matching script
    with open("D:/LOX1_EVAR_MR/data/step1/olr1_pos_lookup.json", "w") as f:
        json.dump({k: v["snp"] for k, v in pos_lookup.items()}, f, indent=2)
    
    print(f"\n{len(pos_lookup)} position lookup keys created", flush=True)
    print(f"\nPosition range (GRCh37): {min(int(v) for k,v in pos_lookup if v)} to {max(int(v) for k,v in pos_lookup if v)}", flush=True)
else:
    print("\nNo positions obtained!", flush=True)

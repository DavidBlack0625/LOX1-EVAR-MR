#!/usr/bin/env python3
"""
Get GRCh37 positions for 31 OLR1 SNPs via Ensembl REST API (batch).
"""
import json, csv, os, urllib.request, ssl, sys, time
sys.stdout.reconfigure(line_buffering=True)

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

# Read SNPs
snps = []
with open("D:/LOX1_EVAR_MR/data/exposure/eqtlgen_OLR1_blood_eqtl.csv", "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        snps.append(row["SNP"])

print(f"Getting GRCh37 positions for {len(snps)} SNPs via Ensembl REST API...", flush=True)

# Ensembl batch POST endpoint
# POST https://rest.ensembl.org/variation/human
# Body: {"ids": ["rs1", "rs2", ...]}
url = "https://rest.ensembl.org/variation/human"
headers = {"Content-Type": "application/json", "Accept": "application/json"}

# Batch in groups of 100 (Ensembl limit)
snp_positions = []
for i in range(0, len(snps), 100):
    batch = snps[i:i+100]
    payload = json.dumps({"ids": batch}).encode("utf-8")
    
    req = urllib.request.Request(url, data=payload, method="POST")
    for k, v in headers.items():
        req.add_header(k, v)
    
    try:
        with urllib.request.urlopen(req, timeout=30, context=ctx) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            
            for snp_id, info in data.items():
                mappings = info.get("mappings", [])
                for m in mappings:
                    if m.get("assembly") == "GRCh37":
                        loc = m.get("location", "")
                        # Format: "12:10310801-10310802"
                        parts = loc.split(":")
                        chrom = parts[0]
                        pos_range = parts[1].split("-")
                        pos = int(pos_range[0])
                        
                        snp_positions.append({
                            "snp": snp_id,
                            "chr_grch37": chrom,
                            "pos_grch37": pos,
                            "allele_string": m.get("allele_string", ""),
                            "assembly": "GRCh37"
                        })
                        print(f"  {snp_id}: chr{chrom}:{pos} ({m.get('allele_string', '')})", flush=True)
                        break  # Take first GRCh37 mapping
                else:
                    # No GRCh37 mapping, try GRCh38
                    for m in mappings:
                        if m.get("assembly") == "GRCh38":
                            loc = m.get("location", "")
                            parts = loc.split(":")
                            chrom = parts[0]
                            pos = int(parts[1].split("-")[0])
                            snp_positions.append({
                                "snp": snp_id,
                                "chr_grch37": "",  # Need liftover
                                "pos_grch37": "",
                                "chr_grch38": chrom,
                                "pos_grch38": pos,
                                "allele_string": m.get("allele_string", ""),
                                "assembly": "GRCh38_only"
                            })
                            print(f"  {snp_id}: chr{chrom}:{pos} (GRCh38 only, {m.get('allele_string', '')})", flush=True)
                            break
    except Exception as e:
        print(f"  Error: {e}", flush=True)
    
    time.sleep(1)  # Ensembl rate limit: 15 requests/second

print(f"\nGot positions for {len(snp_positions)}/{len(snps)} SNPs", flush=True)

# Save to CSV
if snp_positions:
    # Merge with eQTL data
    eqtl_data = {}
    with open("D:/LOX1_EVAR_MR/data/exposure/eqtlgen_OLR1_blood_eqtl.csv", "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            eqtl_data[row["SNP"]] = row
    
    output = []
    for sp in snp_positions:
        eqtl = eqtl_data.get(sp["snp"], {})
        sp["olr1_beta"] = eqtl.get("beta", "")
        sp["olr1_se"] = eqtl.get("se", "")
        sp["olr1_pval"] = eqtl.get("pval", "")
        sp["olr1_ea"] = eqtl.get("effect_allele", "")
        sp["olr1_oa"] = eqtl.get("other_allele", "")
        output.append(sp)
    
    # Write
    fieldnames = ["snp", "chr_grch37", "pos_grch37", "chr_grch38", "pos_grch38", 
                  "allele_string", "assembly", "olr1_beta", "olr1_se", "olr1_pval", "olr1_ea", "olr1_oa"]
    # Remove empty keys
    fieldnames = [f for f in fieldnames if any(f in row for row in output)]
    
    with open("D:/LOX1_EVAR_MR/data/step1/olr1_snp_positions.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(output)
    
    print(f"Saved to olr1_snp_positions.csv", flush=True)
    
    # Print chr12 summary
    chr12_snps = [s for s in snp_positions if s.get("chr_grch37") == "12"]
    print(f"\nChr12 (GRCh37) SNPs: {len(chr12_snps)}", flush=True)
    for s in chr12_snps:
        print(f"  {s['snp']}: 12:{s['pos_grch37']} ({s['allele_string']})", flush=True)

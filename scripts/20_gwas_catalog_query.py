#!/usr/bin/env python3
"""
Step 1: Use GWAS Catalog API (which works!) to find IL6R/CRP associations
for our OLR1 eQTL SNPs. Also search for IL6R protein level studies.
"""
import json, csv, os, time, urllib.request, ssl, sys
sys.stdout.reconfigure(line_buffering=True)

DATA_DIR = "D:/LOX1_EVAR_MR/data"
RESULTS_DIR = "D:/LOX1_EVAR_MR/results"

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

def api_get(url, timeout=20):
    req = urllib.request.Request(url)
    req.add_header("User-Agent", "Mozilla/5.0 Research")
    req.add_header("Accept", "application/json")
    try:
        with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
            raw = resp.read().decode("utf-8")
            try:
                return json.loads(raw)
            except:
                return raw
    except Exception as e:
        return {"error": str(e)}

# Read OLR1 SNPs
snps = []
with open(f"{DATA_DIR}/exposure/eqtlgen_OLR1_blood_eqtl.csv", "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        snps.append(row["SNP"])
print(f"Loaded {len(snps)} OLR1 cis-eQTL SNPs", flush=True)

# ============================================================
# 1. Search GWAS Catalog for IL6R protein studies with different terms
# ============================================================
print("\n" + "="*60, flush=True)
print("1. GWAS Catalog - search for IL6R studies", flush=True)
print("="*60, flush=True)

il6r_studies = []
for search_term in ["interleukin-6", "IL-6", "soluble+IL-6+receptor", "sIL6R", "IL6R"]:
    url = f"https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByDiseaseTrait?diseaseTrait={search_term}"
    result = api_get(url, timeout=15)
    if isinstance(result, dict) and "_embedded" in result:
        studies = result["_embedded"].get("studies", [])
        print(f"  '{search_term}': {len(studies)} studies", flush=True)
        for s in studies:
            acc = s.get("accessionId", "")
            trait = s.get("diseaseTrait", "")
            if isinstance(trait, dict):
                trait = trait.get("trait", "")
            iss = s.get("initialSampleSize", "")
            pub = s.get("publicationInfo", {})
            pub_title = pub.get("title", "") if isinstance(pub, dict) else ""
            pub_year = pub.get("publicationDate", "")[:4] if isinstance(pub, dict) else ""
            il6r_studies.append({
                "acc": acc, "trait": trait, "sample": iss,
                "title": pub_title, "year": pub_year
            })
            print(f"    {acc}: {trait} | {pub_title[:60]} | {iss[:50]} | {pub_year}", flush=True)
    time.sleep(0.3)

# Deduplicate
seen = set()
unique_studies = []
for s in il6r_studies:
    if s["acc"] not in seen:
        seen.add(s["acc"])
        unique_studies.append(s)
print(f"\n  Total unique IL6R studies: {len(unique_studies)}", flush=True)

# ============================================================
# 2. Query each OLR1 SNP in GWAS Catalog for ALL associations
# ============================================================
print("\n" + "="*60, flush=True)
print("2. GWAS Catalog - query each OLR1 SNP for associations", flush=True)
print("="*60, flush=True)

all_snp_assocs = []
protein_hits = []

for i, snp in enumerate(snps):
    url = f"https://www.ebi.ac.uk/gwas/rest/api/associations/search/findByRsId?rsId={snp}"
    result = api_get(url, timeout=15)
    
    if isinstance(result, dict) and "_embedded" in result:
        assoc_list = result["_embedded"].get("associations", [])
        
        for assoc in assoc_list:
            # Extract trait info
            efo_traits = assoc.get("efoTraits", {})
            if isinstance(efo_traits, dict):
                traits_list = efo_traits.get("_embedded", {}).get("efoTraits", [])
                trait_names = [t.get("trait", "") for t in traits_list if isinstance(t, dict)]
            else:
                trait_names = []
            
            # Extract study info
            study_data = assoc.get("study", {})
            if isinstance(study_data, dict):
                study_embedded = study_data.get("_embedded", {}).get("studies", [])
                study_info = study_embedded[0] if study_embedded else {}
            else:
                study_info = {}
            
            reported_trait = study_info.get("diseaseTrait", "")
            if isinstance(reported_trait, dict):
                reported_trait = reported_trait.get("trait", "")
            
            pval = assoc.get("pvalue")
            effect_allele = None
            other_allele = None
            beta = None
            or_val = None
            
            # Extract effect size info
            loci = assoc.get("loci", {})
            if isinstance(loci, dict):
                loci_list = loci.get("_embedded", {}).get("loci", [])
                for locus in loci_list:
                    strongest = locus.get("strongestRiskAlleles", [])
                    for sa in strongest:
                        if isinstance(sa, dict):
                            ra = sa.get("riskAlleleName", "")
                            if ra:
                                # Extract allele from format like "rs123-A"
                                parts = ra.split("-")
                                if len(parts) >= 2:
                                    effect_allele = parts[-1]
            
            snp_assoc = {
                "snp": snp,
                "reported_trait": reported_trait,
                "efo_traits": "; ".join(trait_names),
                "pvalue": pval,
                "accession": study_info.get("accessionId", ""),
                "sample_size": study_info.get("initialSampleSize", ""),
                "effect_allele": effect_allele,
            }
            all_snp_assocs.append(snp_assoc)
            
            # Check if this is a protein-level trait
            trait_combined = (reported_trait + " " + " ".join(trait_names)).lower()
            if any(kw in trait_combined for kw in ["il-6", "il6", "interleukin-6", "crp", "c-reactive", "protein"]):
                protein_hits.append(snp_assoc)
                print(f"  *** {snp}: {reported_trait} (p={pval})", flush=True)
    
    if (i + 1) % 5 == 0:
        print(f"  Processed {i+1}/{len(snps)} SNPs, total associations: {len(all_snp_assocs)}, protein hits: {len(protein_hits)}", flush=True)
    time.sleep(0.3)

print(f"\n  Total associations across all SNPs: {len(all_snp_assocs)}", flush=True)
print(f"  IL6R/CRP/protein hits: {len(protein_hits)}", flush=True)

# Save all associations
if all_snp_assocs:
    with open(f"{DATA_DIR}/step1/gwas_catalog_all_snp_associations.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(all_snp_assocs[0].keys()))
        w.writeheader()
        w.writerows(all_snp_assocs)
    print(f"  Saved all associations to gwas_catalog_all_snp_associations.csv", flush=True)

if protein_hits:
    with open(f"{DATA_DIR}/step1/gwas_catalog_protein_hits.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(protein_hits[0].keys()))
        w.writeheader()
        w.writerows(protein_hits)
    print(f"  Saved protein hits to gwas_catalog_protein_hits.csv", flush=True)

# ============================================================
# 3. Try to get CRP GWAS full summary stats
# ============================================================
print("\n" + "="*60, flush=True)
print("3. Download CRP GWAS summary stats", flush=True)
print("="*60, flush=True)

# GCST000430 is the largest CRP study found (n=17,967)
# Check if full summary stats are available
for acc in ["GCST000430", "GCST000179", "GCST002086"]:
    url = f"https://www.ebi.ac.uk/gwas/rest/api/studies/{acc}"
    result = api_get(url, timeout=15)
    if isinstance(result, dict) and "accessionId" in result:
        # Check for full summary stats URL
        ss_url = None
        if "fullSummaryStats" in result:
            ss_url = result["fullSummaryStats"]
        elif "summaryStatsFile" in result:
            ss_url = result["summaryStatsFile"]
        
        # Check associations count
        assoc_count = result.get("associationCount", 0)
        trait = result.get("diseaseTrait", "")
        if isinstance(trait, dict):
            trait = trait.get("trait", "")
        
        print(f"  {acc}: {trait} | associations={assoc_count} | ss_url={ss_url}", flush=True)
        
        # Try the GWAS Catalog summary stats download path
        # Format: https://www.ebi.ac.uk/gwas/download/summary-statistics/{acc}
        download_url = f"https://www.ebi.ac.uk/gwas/download/summary-statistics/{acc}"
        print(f"    Trying: {download_url}", flush=True)
        dl_result = api_get(download_url, timeout=15)
        if isinstance(dl_result, dict) and "error" in dl_result:
            print(f"    Error: {dl_result['error']}", flush=True)
        elif isinstance(dl_result, str) and "href" in dl_result:
            print(f"    Found download link: {dl_result[:200]}", flush=True)
        elif isinstance(dl_result, dict) and "_links" in dl_result:
            links = dl_result["_links"]
            if "self" in links:
                print(f"    Self link: {links['self'].get('href', '')}", flush=True)
    time.sleep(0.3)

# ============================================================
# 4. Check the newer CRP studies (larger sample sizes)
# ============================================================
print("\n" + "="*60, flush=True)
print("4. Search for larger CRP/IL6R GWAS", flush=True)
print("="*60, flush=True)

# Search for CRP with different terms
for term in ["CRP+level", "C-reactive+protein+level", "high+sensitivity+CRP"]:
    url = f"https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByDiseaseTrait?diseaseTrait={term}"
    result = api_get(url, timeout=15)
    if isinstance(result, dict) and "_embedded" in result:
        studies = result["_embedded"].get("studies", [])
        print(f"  '{term}': {len(studies)} studies", flush=True)
        for s in studies[:5]:
            acc = s.get("accessionId", "")
            trait = s.get("diseaseTrait", "")
            if isinstance(trait, dict):
                trait = trait.get("trait", "")
            iss = s.get("initialSampleSize", "")
            ac = s.get("associationCount", 0)
            print(f"    {acc}: {trait} | n={iss[:50]} | assocs={ac}", flush=True)
    time.sleep(0.3)

# Also search for the specific GCST IDs from our earlier analysis
print("\n  Checking GCST90018950 (CRP, n=436,491 from earlier search)...", flush=True)
url = "https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90018950"
result = api_get(url, timeout=15)
if isinstance(result, dict) and "accessionId" in result:
    trait = result.get("diseaseTrait", "")
    if isinstance(trait, dict):
        trait = trait.get("trait", "")
    print(f"    Found: {trait} | assocs={result.get('associationCount', 0)}", flush=True)
    # Check for summary stats
    links = result.get("_links", {})
    for key, val in links.items():
        if isinstance(val, dict) and "href" in val:
            print(f"    Link {key}: {val['href']}", flush=True)
elif isinstance(result, dict) and "error" in result:
    print(f"    Error: {result['error']}", flush=True)

# ============================================================
# 5. Summary and recommendations
# ============================================================
print("\n" + "="*60, flush=True)
print("SUMMARY", flush=True)
print("="*60, flush=True)
print(f"Total OLR1 SNPs queried: {len(snps)}", flush=True)
print(f"Total GWAS Catalog associations found: {len(all_snp_assocs)}", flush=True)
print(f"IL6R/CRP/protein hits: {len(protein_hits)}", flush=True)
print(f"IL6R studies in GWAS Catalog: {len(unique_studies)}", flush=True)

if not protein_hits:
    print("\nNo direct IL6R/CRP associations found for OLR1 eQTL SNPs.", flush=True)
    print("This is consistent with eQTLGen trans-eQTL null result.", flush=True)
    print("\nInterpretation: OLR1 genetic variants do not directly affect", flush=True)
    print("IL6R or CRP levels. The mediation likely occurs at:", flush=True)
    print("  - Post-transcriptional level (protein processing/cleavage)", flush=True)
    print("  - Through intermediate cell types (macrophage activation)", flush=True)
    print("  - Via signaling pathways (LOX-1 -> oxidative stress -> IL6R)", flush=True)
    print("\nRecommendation: Report Step 1 as null (genetic level),", flush=True)
    print("and frame the study as:", flush=True)
    print("  'OLR1 expression increases AAA risk (MR positive)',")
    print("  'IL6R signaling protects against AAA (Step 2 MR)',")
    print("  'But OLR1 does not genetically regulate IL6R (Step 1 null)',")
    print("  'Suggesting parallel rather than sequential pathways'")

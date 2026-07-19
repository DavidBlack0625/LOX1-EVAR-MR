#!/usr/bin/env python3
"""
Step 1 Mediation MR: Query multiple open APIs for OLR1 SNP -> IL6R/CRP protein associations.
No registration required - uses GWAS Catalog, Open Targets Genetics, and PhenoScanner APIs.
"""

import json
import csv
import os
import time
import math
import urllib.request
import urllib.parse
import urllib.error
import ssl

# === Config ===
DATA_DIR = "D:/LOX1_EVAR_MR/data"
RESULTS_DIR = "D:/LOX1_EVAR_MR/results"
os.makedirs(f"{DATA_DIR}/step1", exist_ok=True)
os.makedirs(f"{RESULTS_DIR}/tables", exist_ok=True)

# Read OLR1 eQTL SNPs
snps = []
with open(f"{DATA_DIR}/exposure/eqtlgen_OLR1_blood_eqtl.csv", "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        snps.append(row["SNP"])

print(f"Loaded {len(snps)} OLR1 cis-eQTL SNPs from eQTLGen")
print(f"SNPs: {snps[:5]}...")

# SSL context that doesn't verify (for some APIs)
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

def api_get(url, headers=None, timeout=30):
    """GET request with error handling"""
    req = urllib.request.Request(url)
    if headers:
        for k, v in headers.items():
            req.add_header(k, v)
    req.add_header("User-Agent", "Mozilla/5.0 Research/1.0")
    try:
        with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
            data = resp.read().decode("utf-8")
            return json.loads(data) if data.startswith("{") or data.startswith("[") else data
    except urllib.error.HTTPError as e:
        return {"error": f"HTTP {e.code}", "msg": str(e)}
    except Exception as e:
        return {"error": str(e)}

def api_post(url, payload, headers=None, timeout=30):
    """POST request with JSON payload"""
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(url, data=data, method="POST")
    req.add_header("Content-Type", "application/json")
    req.add_header("User-Agent", "Mozilla/5.0 Research/1.0")
    if headers:
        for k, v in headers.items():
            req.add_header(k, v)
    try:
        with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except Exception as e:
        return {"error": str(e)}

# ============================================================
# Method 1: GWAS Catalog REST API - search by SNP
# ============================================================
print("\n" + "="*60)
print("Method 1: GWAS Catalog REST API")
print("="*60)

gwas_results = []
for snp in snps:
    url = f"https://www.ebi.ac.uk/gwas/rest/api/associations/search/findByRsId?rsId={snp}"
    result = api_get(url, timeout=15)
    if isinstance(result, dict) and "_embedded" in result:
        assoc_list = result["_embedded"].get("associations", [])
        for assoc in assoc_list:
            trait_name = ""
            if "efoTraits" in assoc:
                traits = assoc["efoTraits"].get("_embedded", {}).get("efoTraits", [])
                trait_name = "; ".join([t.get("trait", "") for t in traits])
            
            study = assoc.get("study", {}).get("_embedded", {}).get("studies", [{}])
            study_info = study[0] if study else {}
            pub_title = study_info.get("publicationInfo", {}).get("title", "") if isinstance(study_info, dict) else ""
            
            # Filter for IL6R, CRP, or inflammatory markers
            trait_lower = trait_name.lower()
            title_lower = pub_title.lower()
            if any(kw in trait_lower or kw in title_lower for kw in 
                   ["il-6", "il6", "interleukin-6", "interleukin 6", "crp", "c-reactive", "inflammation"]):
                pval = assoc.get("pvalue", "NA")
                effect = assoc.get("effectSize", "NA")
                gwas_results.append({
                    "snp": snp,
                    "trait": trait_name,
                    "pvalue": pval,
                    "effect_size": effect,
                    "study_title": pub_title[:100],
                    "source": "GWAS_Catalog"
                })
                print(f"  {snp}: {trait_name} (p={pval})")
    time.sleep(0.3)  # Rate limiting

print(f"\nGWAS Catalog: found {len(gwas_results)} IL6R/CRP/inflammation associations")

# ============================================================
# Method 2: Open Targets Genetics GraphQL API
# ============================================================
print("\n" + "="*60)
print("Method 2: Open Targets Genetics GraphQL API")
print("="*60)

# Query for each SNP: what genes/traits does it associate with?
ot_results = []

for snp in snps[:10]:  # Start with first 10 to test
    # First, get the variant ID in Open Targets format
    # Open Targets uses chr_pos_ref_alt_build format
    # We need to find the rsid -> variant mapping
    
    query = """
    query($rsId: String!) {
        search(queryString: "%s", entityNames: ["variant"]) {
            hits {
                id
                name
                entity
            }
        }
    }
    """ % snp
    
    result = api_post(
        "https://api.genetics.opentargets.org/graphql",
        {"query": query}
    )
    
    if "data" in result and result["data"]:
        hits = result["data"].get("search", {}).get("hits", [])
        if hits:
            variant_id = hits[0]["id"]
            print(f"  {snp} -> OT variant ID: {variant_id}")
            
            # Now query variant associations with IL6R and CRP
            # Use the studies associations endpoint
            query2 = """
            query($variantId: String!) {
                variantInfo(variantId: $variantId) {
                    variantId
                    rsId
                    chromosome
                    position
                }
            }
            """
            result2 = api_post(
                "https://api.genetics.opentargets.org/graphql",
                {"query": query2, "variables": {"variantId": variant_id}}
            )
            if "data" in result2:
                print(f"    Position: {result2['data'].get('variantInfo', {})}")
            
            # Query for gene-burden associations at this variant
            query3 = """
            query($variantId: String!) {
                studiesAndLeadVariantsForVariant(variantId: $variantId) {
                    study {
                        studyId
                        traitReported
                        traitCategory
                    }
                    leadVariant {
                        variantId
                        rsId
                    }
                    pval
                    beta
                    oddsRatio
                }
            }
            """
            result3 = api_post(
                "https://api.genetics.opentargets.org/graphql",
                {"query": query3, "variables": {"variantId": variant_id}}
            )
            
            if "data" in result3 and result3["data"]:
                studies = result3["data"].get("studiesAndLeadVariantsForVariant", [])
                il6_crp_studies = [s for s in studies if 
                    any(kw in (s.get("study", {}).get("traitReported", "") or "").lower() for kw in 
                        ["il-6", "il6", "interleukin", "crp", "c-reactive", "inflammatory"])]
                
                if il6_crp_studies:
                    for s in il6_crp_studies:
                        ot_results.append({
                            "snp": snp,
                            "trait": s["study"]["traitReported"],
                            "category": s["study"].get("traitCategory", ""),
                            "pval": s.get("pval"),
                            "beta": s.get("beta"),
                            "or": s.get("oddsRatio"),
                            "study_id": s["study"]["studyId"],
                            "source": "OpenTargets"
                        })
                        print(f"    FOUND: {s['study']['traitReported']} (p={s.get('pval')}, beta={s.get('beta')})")
        else:
            print(f"  {snp}: no variant found in Open Targets")
    time.sleep(0.5)

print(f"\nOpen Targets: found {len(ot_results)} IL6R/CRP associations")

# ============================================================
# Method 3: PhenoScanner V2 API
# ============================================================
print("\n" + "="*60)
print("Method 3: PhenoScanner V2 API")
print("="*60)

ps_results = []
# PhenoScanner API endpoint
ps_url = "http://www.phenoscanner.medschl.cam.ac.uk/api/"

# Query SNPs in batches
for snp in snps[:15]:  # Try first 15
    params = {
        "apisecret": "phenoscanner",
        "rsid": snp,
        "catalog": "GWAS"
    }
    
    query_string = urllib.parse.urlencode(params)
    full_url = f"{ps_url}?{query_string}"
    
    result = api_get(full_url, timeout=20)
    
    if isinstance(result, dict) and "results" in result:
        for r in result["results"]:
            trait = r.get("PHENO", "")
            if any(kw in trait.lower() for kw in ["il-6", "il6", "interleukin-6", "crp", "c-reactive"]):
                ps_results.append({
                    "snp": snp,
                    "trait": trait,
                    "pvalue": r.get("P", ""),
                    "beta": r.get("BETA", ""),
                    "effect_allele": r.get("EFFECT_ALLELE", ""),
                    "other_allele": r.get("OTHER_ALLELE", ""),
                    "study": r.get("STUDY", ""),
                    "source": "PhenoScanner"
                })
                print(f"  {snp}: {trait} (p={r.get('P')})")
    elif isinstance(result, dict) and "error" in result:
        print(f"  {snp}: {result['error']}")
    
    time.sleep(0.5)

print(f"\nPhenoScanner: found {len(ps_results)} IL6R/CRP associations")

# ============================================================
# Method 4: GWAS Catalog - search for IL6R pQTL studies directly
# ============================================================
print("\n" + "="*60)
print("Method 4: GWAS Catalog - search for IL6R/CRP pQTL studies")
print("="*60)

# Search for studies about IL6R protein levels
for trait_query in ["interleukin-6+receptor", "IL-6+receptor+protein", "C-reactive+protein", "CRP+protein"]:
    url = f"https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByDiseaseTrait?diseaseTrait={trait_query}"
    result = api_get(url, timeout=20)
    
    if isinstance(result, dict) and "_embedded" in result:
        studies = result["_embedded"].get("studies", [])
        for study in studies[:5]:
            acc = study.get("accessionId", "")
            title = study.get("diseaseTrait", "")
            pub = study.get("publicationInfo", {}).get("title", "") if isinstance(study.get("publicationInfo"), dict) else ""
            sample = study.get("initialSampleSize", "")
            print(f"  Study: {acc} | {title} | {pub[:60]} | {sample[:50]}")
    time.sleep(0.3)

# ============================================================
# Method 5: UKB-PPP via GWAS Catalog (Sun et al. 2023 Nature)
# ============================================================
print("\n" + "="*60)
print("Method 5: UKB-PPP - search GWAS Catalog for accession")
print("="*60)

# UKB-PPP study accessions start with GCST90274823+
# Search for "Olink" or "SomaScan" protein studies
url = "https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByDiseaseTrait?diseaseTrait=Olink+protein"
result = api_get(url, timeout=20)
if isinstance(result, dict) and "_embedded" in result:
    studies = result["_embedded"].get("studies", [])
    print(f"Found {len(studies)} Olink protein studies in GWAS Catalog")
    for study in studies[:10]:
        acc = study.get("accessionId", "")
        trait = study.get("diseaseTrait", "")
        print(f"  {acc}: {trait}")
elif isinstance(result, dict) and "error" in result:
    print(f"  Error: {result['error']}")

# Also try searching by publication
url2 = "https://www.ebi.ac.uk/gwas/rest/api/efoTraits/search/findByTraitIri?traitIri=http://www.ebi.ac.uk/efo/EFO_0004760"
result2 = api_get(url2, timeout=20)
print(f"  EFO IL6R trait lookup: {type(result2)}")

# ============================================================
# Method 6: Try IEU OpenGWAS v2 without token (public endpoints)
# ============================================================
print("\n" + "="*60)
print("Method 6: IEU OpenGWAS v2 - public list endpoints")
print("="*60)

# The v2 API might allow listing without token
url = "https://gwas-api.mrcieu.ac.uk/gwasinfo"
result = api_get(url, timeout=20)
if isinstance(result, dict):
    # Search for IL6R/CRP related datasets
    il6_crp_studies = []
    for study_id, info in result.items():
        trait = (info.get("trait", "") or "").lower()
        category = (info.get("category", "") or "").lower()
        if any(kw in trait or kw in category for kw in ["il-6", "il6", "interleukin-6", "crp", "c-reactive", "protein"]):
            il6_crp_studies.append({
                "id": study_id,
                "trait": info.get("trait", ""),
                "category": info.get("category", ""),
                "sample_size": info.get("sample_size", ""),
                "ncase": info.get("ncase", ""),
                "ncontrol": info.get("ncontrol", ""),
            })
    
    print(f"Found {len(il6_crp_studies)} IL6R/CRP/protein datasets in OpenGWAS:")
    for s in il6_crp_studies[:20]:
        print(f"  {s['id']}: {s['trait']} (n={s['sample_size']}, cat={s['category']})")
    
    # Save for reference
    if il6_crp_studies:
        with open(f"{DATA_DIR}/step1/opengwas_il6_crp_studies.csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["id", "trait", "category", "sample_size", "ncase", "ncontrol"])
            writer.writeheader()
            writer.writerows(il6_crp_studies)
        print(f"\n  Saved to {DATA_DIR}/step1/opengwas_il6_crp_studies.csv")
        
        # Try to extract SNPs from the most relevant datasets (without token might fail)
        # But let's try the toptopvariants endpoint
        for study in il6_crp_studies[:3]:
            sid = study["id"]
            if "prot" in sid.lower() or "olink" in sid.lower() or "protein" in sid.lower():
                print(f"\n  Trying to get top variants for {sid}...")
                url2 = f"https://gwas-api.mrcieu.ac.uk/topvariants/{sid}"
                result2 = api_get(url2, timeout=20)
                if isinstance(result2, list) and len(result2) > 0:
                    print(f"    Got {len(result2)} top variants!")
                    # Check if any of our SNPs match
                    for v in result2:
                        if v.get("rsid") in snps:
                            print(f"    *** MATCH: {v['rsid']} in {sid}! beta={v.get('beta')}, p={v.get('p')}")
                elif isinstance(result2, dict) and "error" in result2:
                    print(f"    Error: {result2['error']}")
                time.sleep(0.5)
else:
    print(f"  Error accessing OpenGWAS: {result}")

# ============================================================
# Method 7: Try associating via GWAS Catalog batch query
# ============================================================
print("\n" + "="*60)
print("Method 7: GWAS Catalog batch SNP query for IL6R/CRP")
print("="*60)

# Use the GWAS Catalog associations endpoint with trait filtering
# Search for all associations with IL6R trait
for efo_id in ["EFO_0004760", "EFO_0004474"]:  # IL6R, CRP EFO IDs
    url = f"https://www.ebi.ac.uk/gwas/rest/api/efoTraits/{efo_id}/associations"
    result = api_get(url, timeout=30)
    
    if isinstance(result, dict) and "_embedded" in result:
        assoc_list = result["_embedded"].get("associations", [])
        print(f"  EFO {efo_id}: {len(assoc_list)} total associations")
        
        # Check if any of our SNPs are in the results
        snp_matches = []
        for assoc in assoc_list:
            loci = assoc.get("loci", {}).get("_embedded", {}).get("loci", [])
            for locus in loci:
                strongest = locus.get("strongestRiskAlleles", [])
                for allele in strongest:
                    rsid = allele.get("riskAlleleName", "").replace("*", "")
                    if rsid in snps:
                        pval = assoc.get("pvalue", "NA")
                        snp_matches.append({
                            "snp": rsid,
                            "efo": efo_id,
                            "pvalue": pval
                        })
                        print(f"    *** MATCH: {rsid} (p={pval})")
        
        if snp_matches:
            # Save to file
            with open(f"{DATA_DIR}/step1/gwas_catalog_snp_matches.csv", "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=["snp", "efo", "pvalue", "source"])
                writer.writeheader()
                for m in snp_matches:
                    writer.writerow({**m, "source": "GWAS_Catalog_EFO"})
    elif isinstance(result, dict) and "error" in result:
        print(f"  Error: {result['error']}")
    time.sleep(0.5)

# ============================================================
# Summary
# ============================================================
print("\n" + "="*60)
print("SUMMARY")
print("="*60)
print(f"GWAS Catalog (Method 1): {len(gwas_results)} associations")
print(f"Open Targets (Method 2): {len(ot_results)} associations")
print(f"PhenoScanner (Method 3): {len(ps_results)} associations")
print(f"Total unique: {len(gwas_results) + len(ot_results) + len(ps_results)}")

# Combine all results
all_results = gwas_results + ot_results + ps_results

if all_results:
    with open(f"{RESULTS_DIR}/tables/step1_snp_mediator_lookup.csv", "w", newline="") as f:
        if all_results:
            writer = csv.DictWriter(f, fieldnames=list(all_results[0].keys()))
            writer.writeheader()
            writer.writerows(all_results)
    print(f"\nResults saved to {RESULTS_DIR}/tables/step1_snp_mediator_lookup.csv")
else:
    print("\nNo IL6R/CRP associations found for OLR1 SNPs in any database.")
    print("This supports the hypothesis that OLR1 does not directly regulate IL6R/CRP")
    print("at the genetic level. The mediation may occur through downstream mechanisms.")

print("\n" + "="*60)
print("Next best approach: Use literature-reported IL6R cis-pQTL as instruments")
print("for a REVERSE Mendelian randomization: IL6R -> OLR1 expression")
print("="*60)

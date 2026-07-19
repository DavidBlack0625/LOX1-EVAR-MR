#!/usr/bin/env python3
"""
Step 1: Focused query - find IL6R/CRP pQTL data sources.
Tests APIs one at a time with immediate output.
"""
import json, csv, os, time, urllib.request, urllib.parse, ssl, sys

# Force unbuffered output
sys.stdout.reconfigure(line_buffering=True)

DATA_DIR = "D:/LOX1_EVAR_MR/data"
RESULTS_DIR = "D:/LOX1_EVAR_MR/results"
os.makedirs(f"{DATA_DIR}/step1", exist_ok=True)

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

def api_post(url, payload, timeout=20):
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(url, data=data, method="POST")
    req.add_header("Content-Type", "application/json")
    req.add_header("User-Agent", "Mozilla/5.0 Research")
    try:
        with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except Exception as e:
        return {"error": str(e)}

# ============================================================
# 1. OpenGWAS v2 - list all datasets (no token needed for listing)
# ============================================================
print("="*60, flush=True)
print("1. OpenGWAS v2 - listing datasets...", flush=True)
print("="*60, flush=True)

result = api_get("https://gwas-api.mrcieu.ac.uk/gwasinfo", timeout=30)
if isinstance(result, dict) and "error" not in result:
    print(f"  Total datasets: {len(result)}", flush=True)
    
    # Filter for protein/proteomics datasets
    protein_studies = []
    il6_crp_studies = []
    
    for sid, info in result.items():
        trait = (info.get("trait", "") or "").lower()
        cat = (info.get("category", "") or "").lower()
        
        # Check for IL6R/CRP specifically
        if any(kw in trait for kw in ["il-6", "il6", "interleukin-6", "interleukin 6"]):
            il6_crp_studies.append({"id": sid, "trait": info.get("trait",""), "cat": cat, "n": info.get("sample_size","")})
        if any(kw in trait for kw in ["c-reactive", "crp"]):
            il6_crp_studies.append({"id": sid, "trait": info.get("trait",""), "cat": cat, "n": info.get("sample_size","")})
        
        # Check for protein panels (UKB-PPP, SomaScan, etc.)
        if any(kw in cat or kw in trait for kw in ["protein", "proteomic", "olink", "somascan", "pqtl"]):
            protein_studies.append({"id": sid, "trait": info.get("trait",""), "cat": cat, "n": info.get("sample_size","")})
    
    print(f"  IL6R/CRP specific datasets: {len(il6_crp_studies)}", flush=True)
    for s in il6_crp_studies:
        print(f"    {s['id']}: {s['trait']} (cat={s['cat']}, n={s['n']})", flush=True)
    
    print(f"  Protein panel datasets: {len(protein_studies)}", flush=True)
    for s in protein_studies[:10]:
        print(f"    {s['id']}: {s['trait']} (cat={s['cat']}, n={s['n']})", flush=True)
    
    # Save protein studies for reference
    if protein_studies:
        with open(f"{DATA_DIR}/step1/opengwas_protein_studies.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["id","trait","cat","n"])
            w.writeheader()
            w.writerows(protein_studies)
        print(f"  Saved {len(protein_studies)} protein studies to CSV", flush=True)
    
    # Try to extract IL6R/CRP SNPs from the specific datasets
    # Even without a token, some endpoints may work
    for study in il6_crp_studies:
        sid = study["id"]
        print(f"\n  Trying to extract top variants from {sid}...", flush=True)
        tv_result = api_get(f"https://gwas-api.mrcieu.ac.uk/topvariants/{sid}", timeout=20)
        if isinstance(tv_result, list) and len(tv_result) > 0:
            print(f"    Got {len(tv_result)} top variants!", flush=True)
            # Save them
            with open(f"{DATA_DIR}/step1/opengwas_{sid}_topvariants.csv", "w", newline="") as f:
                if tv_result:
                    w = csv.DictWriter(f, fieldnames=list(tv_result[0].keys()))
                    w.writeheader()
                    w.writerows(tv_result)
            print(f"    Saved to opengwas_{sid}_topvariants.csv", flush=True)
        elif isinstance(tv_result, dict) and "error" in tv_result:
            print(f"    Error: {tv_result['error']}", flush=True)
        else:
            print(f"    Unexpected result type: {type(tv_result)}", flush=True)
        time.sleep(0.5)
else:
    print(f"  Error: {result}", flush=True)

# ============================================================
# 2. GWAS Catalog - search for IL6R/CRP protein level studies
# ============================================================
print("\n" + "="*60, flush=True)
print("2. GWAS Catalog - IL6R/CRP protein studies...", flush=True)
print("="*60, flush=True)

# Search by reported trait
for trait_search in ["interleukin-6+receptor", "C-reactive+protein"]:
    url = f"https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByDiseaseTrait?diseaseTrait={trait_search}"
    result = api_get(url, timeout=20)
    
    if isinstance(result, dict) and "_embedded" in result:
        studies = result["_embedded"].get("studies", [])
        print(f"\n  '{trait_search}': {len(studies)} studies found", flush=True)
        for s in studies[:8]:
            acc = s.get("accessionId", "")
            dt = s.get("diseaseTrait", "")
            iss = s.get("initialSampleSize", "")
            pub = s.get("publicationInfo", {})
            pub_title = pub.get("title", "") if isinstance(pub, dict) else ""
            print(f"    {acc}: {dt} | {pub_title[:60]} | {iss[:50]}", flush=True)
    elif isinstance(result, dict) and "error" in result:
        print(f"  Error for '{trait_search}': {result['error']}", flush=True)
    else:
        print(f"  No results for '{trait_search}'", flush=True)
    time.sleep(0.3)

# ============================================================
# 3. Open Targets Genetics - chr12 OLR1 region -> IL6R associations
# ============================================================
print("\n" + "="*60, flush=True)
print("3. Open Targets Genetics - OLR1 region -> IL6R...", flush=True)
print("="*60, flush=True)

# Query the chr12:10,158,301-10,172,138 region for associations with IL6R (ENSG00000160712)
# and CRP (ENSG00000132693)
query = """
{
  regionInfo(chromosome: "12", start: 10158301, end: 10172138) {
    chromosome
    start
    end
  }
}
"""
result = api_post("https://api.genetics.opentargets.org/graphql", {"query": query}, timeout=20)
if "data" in result:
    print(f"  Region confirmed: {result['data']['regionInfo']}", flush=True)

# Query for genes in this region and their associations
query2 = """
{
  genesInRegion(chromosome: "12", start: 10158301, end: 10172138) {
    gene {
      id
      symbol
      bioType
    }
  }
}
"""
result2 = api_post("https://api.genetics.opentargets.org/graphql", {"query": query2}, timeout=20)
if "data" in result2:
    genes = result2["data"].get("genesInRegion", [])
    for g in genes:
        gene_id = g["gene"]["id"]
        symbol = g["gene"]["symbol"]
        print(f"  Gene in region: {symbol} ({gene_id})", flush=True)

# Query for lead variants in chr12 OLR1 region and their study associations
query3 = """
{
  studiesAndLeadVariantsForGeneByLd(geneId: "ENSG00000173391", ldPopulation: "EUR") {
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
result3 = api_post("https://api.genetics.opentargets.org/graphql", {"query": query3}, timeout=30)
if "data" in result3 and result3["data"]:
    studies = result3["data"].get("studiesAndLeadVariantsForGeneByLd", [])
    print(f"\n  Total study associations for OLR1 gene: {len(studies)}", flush=True)
    
    # Filter for IL6R/CRP/inflammation
    il6_crp = [s for s in studies if 
        any(kw in (s["study"]["traitReported"] or "").lower() for kw in 
            ["il-6", "il6", "interleukin-6", "crp", "c-reactive", "inflammatory", "inflammation"])]
    
    print(f"  IL6R/CRP/inflammation associations: {len(il6_crp)}", flush=True)
    for s in il6_crp:
        print(f"    {s['study']['traitReported']} | rsID={s['leadVariant']['rsId']} | p={s['pval']} | beta={s.get('beta')} | OR={s.get('oddsRatio')}", flush=True)
    
    # Also check for any protein-level associations
    protein_assoc = [s for s in studies if
        any(kw in (s["study"]["traitReported"] or "").lower() for kw in
            ["protein", "olink", "somascan", "plasma"])]
    print(f"  Protein-level associations: {len(protein_assoc)}", flush=True)
    for s in protein_assoc[:10]:
        print(f"    {s['study']['traitReported']} | rsID={s['leadVariant']['rsId']} | p={s['pval']} | beta={s.get('beta')}", flush=True)
else:
    print(f"  Error or no data: {result3}", flush=True)

# ============================================================
# 4. Directly check if rs2228145 (IL6R functional SNP) affects OLR1
# ============================================================
print("\n" + "="*60, flush=True)
print("4. Reverse check: IL6R SNP rs2228145 -> OLR1 expression", flush=True)
print("="*60, flush=True)

# Query Open Targets for rs2228145 associations
query4 = """
{
  search(queryString: "rs2228145", entityNames: ["variant"]) {
    hits {
      id
      name
      entity
    }
  }
}
"""
result4 = api_post("https://api.genetics.opentargets.org/graphql", {"query": query4}, timeout=20)
if "data" in result4 and result4["data"]:
    hits = result4["data"]["search"]["hits"]
    if hits:
        variant_id = hits[0]["id"]
        print(f"  rs2228145 -> OT variant ID: {variant_id}", flush=True)
        
        # Get all study associations for this variant
        query5 = """
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
        result5 = api_post("https://api.genetics.opentargets.org/graphql", 
                          {"query": query5, "variables": {"variantId": variant_id}}, timeout=30)
        
        if "data" in result5 and result5["data"]:
            studies = result5["data"].get("studiesAndLeadVariantsForVariant", [])
            print(f"  Total associations for rs2228145: {len(studies)}", flush=True)
            
            # Look for OLR1, LOX1, or any chr12 gene associations
            olr1_assoc = [s for s in studies if
                any(kw in (s["study"]["traitReported"] or "").lower() for kw in
                    ["olr1", "lox1", "lox-1", "lectin"])]
            print(f"  OLR1/LOX1 associations: {len(olr1_assoc)}", flush=True)
            
            # Also show top associations for context
            sorted_studies = sorted(studies, key=lambda x: x.get("pval") or 1)
            print(f"  Top 15 associations for rs2228145:", flush=True)
            for s in sorted_studies[:15]:
                print(f"    {s['study']['traitReported']} | p={s.get('pval')} | beta={s.get('beta')} | OR={s.get('oddsRatio')} | cat={s['study'].get('traitCategory','')}", flush=True)
    else:
        print("  rs2228145 not found in Open Targets", flush=True)
else:
    print(f"  Error: {result4}", flush=True)

print("\n" + "="*60, flush=True)
print("DONE - Check results above", flush=True)
print("="*60, flush=True)

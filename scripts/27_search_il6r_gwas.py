"""
Check OpenGWAS API and search for IL6R/CRP pQTL datasets.
Also search GWAS Catalog for IL6R protein GWAS studies.
"""
import json, urllib.request, ssl

# Bypass SSL verification
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

def fetch_json(url, timeout=15):
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
        return json.loads(resp.read().decode('utf-8'))

# 1. Check OpenGWAS API
print("=" * 60)
print("1. OpenGWAS API Status Check")
print("=" * 60)
try:
    d = fetch_json("https://gwas-api.mrcieu.ac.uk/gwasinfo", timeout=20)
    print(f"  STATUS: WORKING ({len(d)} datasets available)")
    
    il6r_found = []
    crp_found = []
    
    for k, v in d.items():
        trait = v.get('trait', '').lower()
        cat = v.get('category', '').lower()
        
        if any(x in trait for x in ['il6r', 'interleukin-6 receptor', 'soluble il-6', 'sil6r']):
            il6r_found.append((k, v.get('trait', ''), v.get('sample_size', '?'), v.get('population', '?')))
        
        if any(x in trait for x in ['c-reactive protein', 'crp ']) and 'protein' in trait:
            crp_found.append((k, v.get('trait', ''), v.get('sample_size', '?'), v.get('population', '?')))
    
    print(f"\n  IL6R-related datasets ({len(il6r_found)}):")
    for id, trait, n, pop in il6r_found[:10]:
        print(f"    {id}: {trait} (n={n}, pop={pop})")
    
    print(f"\n  CRP protein datasets ({len(crp_found)}):")
    for id, trait, n, pop in crp_found[:10]:
        print(f"    {id}: {trait} (n={n}, pop={pop})")
        
except Exception as e:
    print(f"  STATUS: FAILED ({e})")

# 2. Search GWAS Catalog for IL6R protein GWAS
print("\n" + "=" * 60)
print("2. GWAS Catalog: Search for IL6R protein GWAS")
print("=" * 60)
try:
    # Search for IL6R traits
    url = "https://www.ebi.ac.uk/gwas/rest/api/efoTraits/search/findAll?trait=interleukin-6+receptor"
    d = fetch_json(url, timeout=15)
    traits = d.get('_embedded', {}).get('efoTraits', [])
    print(f"  Found {len(traits)} EFO traits for 'interleukin-6 receptor'")
    
    for t in traits[:5]:
        trait_id = t.get('shortForm', '')
        trait_name = t.get('trait', '')
        print(f"    EFO: {trait_id} = {trait_name}")
        
        # Get studies for this trait
        try:
            study_url = f"https://www.ebi.ac.uk/gwas/rest/api/efoTraits/{trait_id}/studies"
            sd = fetch_json(study_url, timeout=15)
            studies = sd.get('_embedded', {}).get('studies', [])
            print(f"      Studies: {len(studies)}")
            for s in studies[:5]:
                acc = s.get('accessionId', '')
                trait_reported = s.get('diseaseTrait', {}).get('trait', '') if s.get('diseaseTrait') else ''
                n = s.get('initialSampleSize', '')
                sumstats = s.get('fullPvalueSet', False)
                print(f"        {acc}: {trait_reported} | n={n} | sumstats={sumstats}")
        except Exception as e:
            print(f"      Error fetching studies: {e}")
            
except Exception as e:
    print(f"  Error: {e}")

# 3. Also search for "soluble IL6 receptor" and "sIL6R"
print("\n--- Also searching 'soluble interleukin-6 receptor' ---")
try:
    url = "https://www.ebi.ac.uk/gwas/rest/api/efoTraits/search/findAll?trait=soluble+interleukin"
    d = fetch_json(url, timeout=15)
    traits = d.get('_embedded', {}).get('efoTraits', [])
    print(f"  Found {len(traits)} EFO traits")
    for t in traits[:5]:
        print(f"    {t.get('shortForm', '')}: {t.get('trait', '')}")
except Exception as e:
    print(f"  Error: {e}")

# 4. Search GWAS Catalog studies directly for IL6R
print("\n--- Direct study search for 'IL6R' or 'interleukin 6 receptor' ---")
try:
    url = "https://www.ebi.ac.uk/gwas/rest/api/studies/search/findAll?trait=interleukin+6+receptor&size=20"
    d = fetch_json(url, timeout=15)
    studies = d.get('_embedded', {}).get('studies', [])
    print(f"  Found {len(studies)} studies")
    for s in studies[:10]:
        acc = s.get('accessionId', '')
        trait_reported = s.get('diseaseTrait', {}).get('trait', '') if s.get('diseaseTrait') else ''
        n = s.get('initialSampleSize', '')
        sumstats = s.get('fullPvalueSet', False)
        print(f"    {acc}: {trait_reported} | n={n} | sumstats={sumstats}")
except Exception as e:
    print(f"  Error: {e}")

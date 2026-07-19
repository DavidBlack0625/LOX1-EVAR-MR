"""
Prepare colocalization analysis data:
1. Re-fetch GTEx eQTL with MAF
2. Process eQTLGen data (compute beta/SE from Z-score)
3. Extract FinnGen regional data for OLR1 ± 1Mb
4. Match SNPs between eQTL and GWAS
5. Save formatted coloc input
"""
import json
import csv
import gzip
import math
import os
import urllib.request

DATA_DIR = "D:/LOX1_EVAR_MR/data"
EXP_DIR = os.path.join(DATA_DIR, "exposure")
OUT_DIR = os.path.join(DATA_DIR, "outcome")
COLOC_DIR = os.path.join(DATA_DIR, "coloc")
os.makedirs(COLOC_DIR, exist_ok=True)

# OLR1 gene location (GRCh38): chr12:10,158,301-10,172,138
# Use ±1Mb window: chr12:9,158,301 - 11,172,138
CHROM = "12"
POS_START = 9158301
POS_END = 111172138  # Actually use ±1Mb from gene center
GENE_CENTER = 10165220  # midpoint of OLR1
WINDOW = 1000000  # ±1Mb
POS_START = GENE_CENTER - WINDOW  # 9165220
POS_END = GENE_CENTER + WINDOW    # 11165220


def fetch_gtex_eqtl_with_maf(tissue):
    """Fetch GTEx eQTL data including MAF from API."""
    url = (
        f"https://gtexportal.org/api/v2/association/singleTissueEqtl"
        f"?gencodeId=ENSG00000173391.8&datasetId=gtex_v8"
        f"&tissueSiteDetailId={tissue}&itemsPerPage=500"
    )
    with urllib.request.urlopen(url, timeout=30) as resp:
        data = json.loads(resp.read().decode())
    return data.get("data", [])


def compute_se_from_pval(beta, pval):
    """Compute SE from beta and p-value using Z-score."""
    if pval <= 0 or pval >= 1:
        return None
    # Two-tailed p-value: |Z| = inv_cdf(1 - p/2)
    from statistics import NormalDist
    z = abs(NormalDist().inv_cdf(1 - pval / 2))
    if z == 0:
        return None
    return abs(beta) / z


def prepare_gtex_data():
    """Prepare GTEx eQTL data with MAF for coloc."""
    results = {}
    for tissue, label, n_samples in [
        ("Artery_Tibial", "Artery_Tibial", 387),
        ("Artery_Coronary", "Artery_Coronary", 118),
    ]:
        records = fetch_gtex_eqtl_with_maf(tissue)
        formatted = []
        for rec in records:
            nes = rec.get("nes", 0)
            pval = rec.get("pValue", 1)
            se = compute_se_from_pval(nes, pval)
            if se is None or se == 0:
                continue

            # Parse variant ID for alleles
            # Format: chr12_10131353_A_C_b38
            vid = rec.get("variantId", "")
            parts = vid.split("_")
            if len(parts) < 5:
                continue
            ref_allele = parts[2]
            alt_allele = parts[3]

            # MAF from API (if available)
            maf = rec.get("maf")
            if maf is None:
                maf = 0.1  # default if not available

            formatted.append({
                "snp": rec.get("snpId", ""),
                "pos": rec.get("pos", 0),
                "beta": nes,
                "se": se,
                "varbeta": se ** 2,
                "maf": maf,
                "N": n_samples,
                "type": "quant",
                "effect_allele": alt_allele,  # NES is for alternate allele
                "other_allele": ref_allele,
                "tissue": tissue,
            })
        results[tissue] = formatted
        print(f"  {tissue}: {len(formatted)} SNPs with MAF and SE")
    return results


def prepare_eqtlgen_data():
    """Prepare eQTLGen eQTL data for coloc."""
    # Read the original eQTLGen significant eQTLs
    eqtlgen_file = os.path.join(EXP_DIR, "eqtlgen_OLR1_blood_eqtl.csv")
    formatted = []

    with open(eqtlgen_file, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            snp = row.get("SNP", "")
            beta_val = float(row.get("beta", 0))
            pval = float(row.get("pval", 1))

            # The beta in our file is actually the Z-score from eQTLGen
            # For standardized expression: beta = Z / sqrt(N), varbeta = 1/N
            # N for eQTLGen = 31684
            N = 31684
            z = beta_val  # This is actually the Z-score
            se = 1.0 / math.sqrt(N)
            beta_actual = z * se
            varbeta = se ** 2

            # MAF - not available in eQTLGen significant file, use 0.1 as default
            # (coloc can work without MAF if we provide beta and varbeta)
            maf = 0.1

            formatted.append({
                "snp": snp,
                "pos": 0,  # not available in our file
                "beta": beta_actual,
                "se": se,
                "varbeta": varbeta,
                "maf": maf,
                "N": N,
                "type": "quant",
                "effect_allele": row.get("effect_allele", ""),
                "other_allele": row.get("other_allele", ""),
                "tissue": "Whole_Blood",
            })

    print(f"  eQTLGen Whole_Blood: {len(formatted)} SNPs")
    return formatted


def extract_finngen_regional():
    """Extract FinnGen SNPs in the OLR1 region (±1Mb)."""
    finngen_gz = "D:/workbuddy工作空间/LOX-1_EVAR_MR/data/outcome/finngen_R12_I9_ABAORTANEUR.gz"

    results = []
    print(f"  Scanning FinnGen for chr12:{POS_START}-{POS_END}...")

    with gzip.open(finngen_gz, "rt") as f:
        header = f.readline().strip().split("\t")
        # Columns: #chrom, pos, ref, alt, rsids, pval, beta, sebeta, af, ...
        col_map = {name: i for i, name in enumerate(header)}

        for line in f:
            parts = line.strip().split("\t")
            chrom = parts[0]
            if chrom != CHROM:
                continue
            pos = int(parts[1])
            if pos < POS_START or pos > POS_END:
                continue

            rsids = parts[4]  # may be semicolon-separated
            pval = float(parts[col_map.get("pval", 6)])
            beta = float(parts[col_map.get("beta", 8)])
            sebeta = float(parts[col_map.get("sebeta", 9)])
            af_alt = float(parts[col_map.get("af_alt", 10)])

            # Use first rsID
            first_rsid = rsids.split(";")[0].strip()
            if not first_rsid:
                continue

            # MAF = min(af_alt, 1 - af_alt)
            maf_val = min(af_alt, 1 - af_alt)
            if maf_val <= 0 or maf_val >= 1:
                continue  # skip monomorphic or invalid SNPs

            results.append({
                "snp": first_rsid,
                "pos": pos,
                "beta": beta,
                "se": sebeta,
                "varbeta": sebeta ** 2,
                "maf": maf_val,
                "pval": pval,
                "af": af_alt,
                "ref": parts[2],
                "alt": parts[3],
                "N": 467545,  # 4439 cases + 463106 controls
                "ncase": 4439,
                "ncontrol": 463106,
                "type": "cc",
                "s": 4439 / 467545,
            })

    print(f"  Found {len(results)} FinnGen SNPs in region")
    return results


def match_snps(eqtl_data, gwas_data):
    """Match SNPs between eQTL and GWAS datasets.
    Use FinnGen MAF for both datasets (MAF should be similar across EUR datasets).
    """
    gwas_by_snp = {d["snp"]: d for d in gwas_data}

    matched = []
    for eq in eqtl_data:
        snp = eq["snp"]
        if snp in gwas_by_snp:
            gw = gwas_by_snp[snp]
            # Use FinnGen MAF for both datasets (more reliable than default 0.1)
            shared_maf = gw["maf"]
            if shared_maf <= 0 or shared_maf >= 1:
                continue  # skip invalid MAF

            matched.append({
                "snp": snp,
                "eqtl_beta": eq["beta"],
                "eqtl_varbeta": eq["varbeta"],
                "eqtl_maf": shared_maf,  # use FinnGen MAF
                "eqtl_N": eq["N"],
                "eqtl_type": "quant",
                "gwas_beta": gw["beta"],
                "gwas_varbeta": gw["varbeta"],
                "gwas_maf": shared_maf,
                "gwas_N": gw["N"],
                "gwas_type": "cc",
                "gwas_s": gw["s"],
                "gwas_pval": gw["pval"],
                "pos": gw.get("pos", 0),
            })
    return matched


def save_coloc_input(matched, tissue_name):
    """Save matched data as CSV for R coloc analysis."""
    if not matched:
        print(f"  WARNING: No matched SNPs for {tissue_name}")
        return None

    filepath = os.path.join(COLOC_DIR, f"coloc_input_{tissue_name}.csv")
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=matched[0].keys())
        writer.writeheader()
        writer.writerows(matched)
    print(f"  Saved {len(matched)} matched SNPs to {filepath}")
    return filepath


def main():
    print("=" * 60)
    print("Preparing colocalization analysis data")
    print("=" * 60)

    # 1. Prepare eQTL data
    print("\n[1] Fetching GTEx eQTL data with MAF...")
    gtex_data = prepare_gtex_data()

    print("\n[2] Preparing eQTLGen data...")
    eqtlgen_data = prepare_eqtlgen_data()

    # 2. Extract FinnGen regional data
    print("\n[3] Extracting FinnGen regional data (OLR1 ± 1Mb)...")
    gwas_data = extract_finngen_regional()

    # Save FinnGen regional data
    if gwas_data:
        fg_file = os.path.join(COLOC_DIR, "finngen_OLR1_region.csv")
        with open(fg_file, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=gwas_data[0].keys())
            writer.writeheader()
            writer.writerows(gwas_data)
        print(f"  Saved to {fg_file}")

    # 3. Match SNPs for each eQTL source
    print("\n[4] Matching SNPs between eQTL and GWAS...")

    all_matched = {}

    for tissue, data in gtex_data.items():
        matched = match_snps(data, gwas_data)
        all_matched[tissue] = matched
        save_coloc_input(matched, tissue)

    # eQTLGen
    matched_eqtlgen = match_snps(eqtlgen_data, gwas_data)
    all_matched["eQTLGen_Blood"] = matched_eqtlgen
    save_coloc_input(matched_eqtlgen, "eQTLGen_Blood")

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for tissue, matched in all_matched.items():
        print(f"  {tissue}: {len(matched)} matched SNPs")

    # Also save all FinnGen regional SNPs for coloc with full region
    # (using only eQTL significant SNPs limits power, but we can also
    # try a broader approach where we use all FinnGen SNPs as the GWAS
    # and only the eQTL SNPs we have)
    print("\nDone! Data ready for coloc.abf analysis.")


if __name__ == "__main__":
    main()

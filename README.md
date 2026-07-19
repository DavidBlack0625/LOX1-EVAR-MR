# OLR1-EVAR-MR: Mendelian Randomization Analysis of OLR1 and Abdominal Aortic Aneurysm

## Overview

This repository contains all analysis code for the Mendelian randomization (MR) study investigating the causal relationship between OLR1 (LOX-1) expression and abdominal aortic aneurysm (AAA) risk, including colocalization and mediation analysis.

## Study Design

- **Exposure**: OLR1 (LOX-1) cis-eQTLs from eQTLGen (whole blood, n=31,684) and GTEx v8 (Artery_Tibial, Artery_Coronary)
- **Outcome**: Abdominal aortic aneurysm (FinnGen R12, I9_ABAORTANEUR, 4,439 cases / 463,106 controls)
- **Mediators**: C-reactive protein (CRP) and interleukin-6 receptor (IL6R)
- **Methods**: Two-sample MR (IVW, MR-Egger, weighted median), coloc.abf, mediation MR (two-step)

## Key Findings

| Analysis | Result |
|----------|--------|
| OLR1 → AAA (total effect) | IVW OR=1.0023, p=0.022 |
| Colocalization | PP.H4 < 2% (negative) |
| CRP → AAA (Step 2) | OR=0.968, p=0.580 (null) |
| IL6R → AAA (Step 2) | OR=0.722, p=2.96×10⁻²⁶ |
| Step 1 (OLR1 → IL6R/CRP) | 0/31 SNP overlap (null) |

## Repository Structure

```
LOX1-EVAR-MR/
├── scripts/                    # Analysis scripts
│   ├── 03_run_mr.R            # GTEx MR analysis (R)
│   ├── 04_mr_eqtlgen.R        # eQTLGen MR analysis (R)
│   ├── 05_prepare_coloc_data.py  # Colocalization data preparation
│   ├── 06_coloc_analysis.R    # Colocalization analysis (R)
│   ├── 12_download_and_extract.py  # FinnGen download + SNP extraction
│   ├── 13_mediation_mr_final.py   # Step 2 mediation MR
│   ├── 14_generate_report.py  # HTML report generation
│   ├── 17_download_trans_eqtl.py  # Step 1 eQTLGen trans-eQTL
│   ├── 26_fix_pvalue_precision.py # P-value precision fix
│   ├── 29_final_data_audit.py # Data completeness audit
│   ├── 33_convert_figures.py  # High-DPI figure conversion
│   └── 34_render_mediation_forest.py  # Mediation forest plot
├── data/                       # Data files (not included - see sources)
│   ├── exposure/              # eQTL exposure data
│   ├── outcome/               # FinnGen AAA GWAS
│   └── mediators/             # CRP/IL6R instruments
├── results/
│   ├── tables/                # Result CSV files
│   └── figures_highres/       # 300 DPI figures (TIFF/PDF/PNG)
└── README.md
```

## Data Sources

| Data | Source | Access |
|------|--------|--------|
| OLR1 eQTL (blood) | eQTLGen (n=31,684) | https://www.eqtlgen.org/ |
| OLR1 eQTL (artery) | GTEx v8 | https://gtexportal.org/ |
| AAA GWAS | FinnGen R12 | https://r12.finngen.fi/ |
| CRP pQTL | Literature (WGHS) | See references in manuscript |
| IL6R pQTL | Literature (IL6R Mendelian randomization) | See references in manuscript |
| trans-eQTL | eQTLGen trans-eQTL | https://www.eqtlgen.org/trans-eqtls.html |

## Software Requirements

- Python 3.13+ (statistics, PIL/Pillow)
- R 4.3+ (TwoSampleMR, coloc, data.table, ggplot2)
- No special hardware required

## Reproducing Results

1. Download FinnGen R12 AAA GWAS from `https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_I9_ABAORTANEUR.gz`
2. Download eQTLGen OLR1 cis-eQTL from `https://www.eqtlgen.org/cis-eqtls.html`
3. Run scripts in numerical order (03 → 04 → 05 → 06 → ... → 34)

## Citation

If you use this code, please cite:
Wei W, Holt PJS, Chen YK. Genetically Predicted OLR1 Expression and Abdominal Aortic Aneurysm Risk: A Mendelian Randomization Study with Colocalization and Mediation Analysis. [Journal name], 2026. (Manuscript in preparation)

## License

MIT License

## Contact

Wei Wei (魏巍) - Department of General Surgery, Chengdu Fifth People's Hospital, China
Email: DavidBlack0625@outlook.com

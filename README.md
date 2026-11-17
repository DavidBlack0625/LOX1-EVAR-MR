# OLR1-EVAR-MR: Mendelian Randomization Analysis of OLR1 and Abdominal Aortic Aneurysm

## Overview

This repository contains all analysis code for the Mendelian randomization (MR) study examining the genetic association between OLR1 (LOX-1) expression and abdominal aortic aneurysm (AAA) risk, including colocalization and a parallel-pathway analysis of CRP and IL6R.

## Study Design

- **Exposure**: OLR1 (LOX-1) cis-eQTLs from eQTLGen (whole blood, n=31,684) and GTEx v8 (Artery_Tibial, Artery_Coronary)
- **Outcome**: Abdominal aortic aneurysm (FinnGen R12, I9_ABAORTANEUR, 4,439 cases / 463,106 controls)
- **Parallel pathways**: C-reactive protein (CRP) and interleukin-6 receptor (IL6R), analysed as independent exposures rather than as mediators
- **Methods**: Two-sample MR (IVW, MR-Egger, weighted median, simple mode, weighted mode), coloc.abf, two-step trans-eQTL screen to test genetic overlap between OLR1 instruments and IL6R/CRP

## Instrument selection

OLR1 cis-eQTLs (within ±1 Mb of ENSG00000173391) were carried forward **without a p-value filter** and clumped for independence (`clump_kb = 10000`, `clump_r2 = 0.001`, 1000 Genomes Phase 3 EUR reference), yielding 31 variants; 24 remained after harmonization with the outcome. Nine of the 31 reached conventional genome-wide significance (p < 5×10⁻⁸); all were retained because per-variant effects on OLR1 expression are modest and restricting to genome-wide significant variants would have substantially reduced instrument strength.

## Effect-size scale

eQTLGen distributes **Z-scores**, not beta coefficients with standard errors. The Z-score was therefore used as the effect estimate with the standard error set to 1, a convention in eQTL-based MR; primary estimates are expressed **per 1-unit increase in the OLR1 eQTL Z-score**. To aid interpretation, per-standard-deviation effect sizes were approximated as

$$\beta_{\text{per SD}} = \frac{Z}{\sqrt{2p(1-p)N}}, \quad N = 31{,}684$$

using 1000 Genomes allele frequencies; 23 of the 24 harmonised variants could be rescaled (rs2896042 is tri-allelic and lacks allele-frequency data).

## Key Findings

| Analysis | Result |
|----------|--------|
| OLR1 → AAA | IVW OR = 1.0023 per 1-unit eQTL Z-score (95% CI 1.0003–1.0043), p = 0.022; equivalent to **OR 1.29 per SD (95% CI 1.03–1.62), p = 0.028** |
| Colocalization | PP.H4 < 2% (negative; PP.H1 > 75%) across all three eQTL sources |
| CRP → AAA | OR = 0.968 (95% CI 0.864–1.085), p = 0.580 (null) |
| IL6R → AAA | OR = 0.722 per unit genetically predicted sIL6R (95% CI 0.680–0.767), p = 2.96×10⁻²⁶; per-Ala358-allele OR 0.86 (95% CI 0.82–0.90) |
| Trans-eQTL screen (OLR1 → IL6R/CRP) | 0/31 SNP overlap (null; the two pathways are genetically independent) |

## Repository Structure

```
LOX1-EVAR-MR/
├── scripts/                          # Analysis scripts
│   ├── 03_run_mr.R                   # GTEx MR analysis (R)
│   ├── 04_mr_eqtlgen.R               # eQTLGen MR analysis (R)
│   ├── 05_prepare_coloc_data.py      # Colocalization data preparation
│   ├── 06_coloc_analysis.R           # Colocalization analysis (R)
│   ├── 12_download_and_extract.py    # FinnGen download + SNP extraction
│   ├── 13_mediation_mr_final.py      # Step 2 parallel-pathway MR
│   ├── 14_generate_report.py         # HTML report generation
│   ├── 17_download_trans_eqtl.py     # Step 1 eQTLGen trans-eQTL
│   ├── 26_fix_pvalue_precision.py    # P-value precision fix
│   ├── 29_final_data_audit.py        # Data completeness audit
│   ├── 33_convert_figures.py         # High-DPI figure conversion
│   └── 34_render_mediation_forest.py # Parallel-pathway forest plot
├── data/                             # Data files (not included - see sources)
│   ├── exposure/                     # eQTL exposure data
│   ├── outcome/                      # FinnGen AAA GWAS
│   └── mediators/                    # CRP/IL6R instruments
├── results/
│   ├── tables/                       # Result CSV files
│   └── figures_highres/              # 300 DPI figures (TIFF/PDF/PNG)
└── README.md
```

> File names such as `13_mediation_mr_final.py` and `34_render_mediation_forest.py` predate the change in framing and are kept for reproducibility; the analysis they perform is a two-step parallel-pathway comparison, not a mediation decomposition.

## Data Sources

| Data | Source | Access |
|------|--------|--------|
| OLR1 eQTL (blood) | eQTLGen (n=31,684) | https://www.eqtlgen.org/ |
| OLR1 eQTL (artery) | GTEx v8 | https://gtexportal.org/ |
| AAA GWAS | FinnGen R12 | https://r12.finngen.fi/ |
| CRP variants | Wensley F, et al. BMJ. 2011;342:d548 (4 CRP variants) | See reference 13 in manuscript |
| IL6R variant | rs2228145 (Asp358Ala); Harrison SC, et al. Eur Heart J. 2013;34(48):3707-3716 | See references 22–24 in manuscript |
| trans-eQTL | eQTLGen trans-eQTL | https://www.eqtlgen.org/trans-eqtls.html |

## Software Requirements

- Python 3.10+ (pandas, numpy, scipy, Pillow)
- R 4.3+ (TwoSampleMR, coloc, data.table, ggplot2)
- No special hardware required

## Reproducing Results

1. Download the FinnGen R12 AAA GWAS from `https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_I9_ABAORTANEUR.gz`
2. Download eQTLGen OLR1 cis-eQTLs from `https://www.eqtlgen.org/cis-eqtls.html`
3. Run scripts in numerical order (03 → 04 → 05 → 06 → … → 34)

## Known issues

- `12_download_and_extract.py` selects the lead SNP with `max(..., key=pval)`, which returns the *least* significant variant; this should be `min`. The OLR1 primary analysis (script 04) is unaffected.
- Results in `results/tables/` were produced before the per-SD rescaling described above; the Z-score scale values are those in the manuscript's Table I.

## Citation

If you use this code, please cite:

Wei W. Genetically Predicted OLR1 Expression and Abdominal Aortic Aneurysm Risk: A Mendelian Randomization Study with Colocalization and Parallel-Pathway Analysis. *Journal of Vascular Surgery: Vascular Science*, 2026. (Manuscript under review)

## License

MIT License

## Contact

Wei Wei (魏巍) — Department of General Surgery, Chengdu Fifth People's Hospital, China
Email: DavidBlack0625@outlook.com

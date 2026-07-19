const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  AlignmentType, BorderStyle, WidthType, ShadingType,
  PageNumber, PageBreak, HeadingLevel, Header, Footer
} = require("docx");

const A4 = { width: 11906, height: 16838 };
const MARGINS = { top: 1440, right: 1080, bottom: 1440, left: 1080 };
const CONTENT_WIDTH = A4.width - MARGINS.left - MARGINS.right;

const border = { style: BorderStyle.SINGLE, size: 1, color: "AAAAAA" };
const borders = { top: border, bottom: border, left: border, right: border };

function cell(text, opts = {}) {
  const isHeader = opts.header || false;
  return new TableCell({
    borders,
    width: { size: opts.width || 1500, type: WidthType.DXA },
    shading: isHeader ? { fill: "D9E2F3", type: ShadingType.CLEAR } : undefined,
    margins: { top: 30, bottom: 30, left: 60, right: 60 },
    children: [new Paragraph({
      children: [new TextRun({ text: String(text), size: 18, bold: isHeader, font: "Times New Roman" })],
      alignment: opts.align || AlignmentType.CENTER,
      spacing: { after: 0, line: 240 },
    })],
  });
}

function makeTable(headers, rows, colWidths) {
  const widths = colWidths || headers.map(() => Math.floor(CONTENT_WIDTH / headers.length));
  const headerRow = new TableRow({
    tableHeader: true,
    children: headers.map((h, i) => cell(h, { header: true, width: widths[i] })),
  });
  const dataRows = rows.map(row =>
    new TableRow({ children: row.map((c, i) => cell(c, { width: widths[i] })) })
  );
  return new Table({
    width: { size: CONTENT_WIDTH, type: WidthType.DXA },
    columnWidths: widths,
    rows: [headerRow, ...dataRows],
  });
}

function heading(text, level) {
  const sizes = { 1: 26, 2: 24 };
  const levels = { 1: HeadingLevel.HEADING_1, 2: HeadingLevel.HEADING_2 };
  return new Paragraph({
    heading: levels[level],
    spacing: { before: 240, after: 120 },
    children: [new TextRun({ text, bold: true, size: sizes[level], font: "Times New Roman" })],
  });
}

function body(text) {
  return new Paragraph({
    children: [new TextRun({ text, size: 22, font: "Times New Roman" })],
    spacing: { after: 80, line: 320 },
    alignment: AlignmentType.JUSTIFIED,
  });
}

function tableTitle(text) {
  return new Paragraph({
    spacing: { before: 200, after: 60 },
    children: [new TextRun({ text, bold: true, size: 20, font: "Times New Roman" })],
  });
}

function spacer() {
  return new Paragraph({ spacing: { after: 80 }, children: [] });
}

// Read CSV data
const singleSnpCsv = fs.readFileSync("D:/LOX1_EVAR_MR/results/tables/single_snp_eQTLGen_Blood.csv", "utf-8");
const looCsv = fs.readFileSync("D:/LOX1_EVAR_MR/results/tables/leaveoneout_eQTLGen_Blood.csv", "utf-8");
const medOutcomesCsv = fs.readFileSync("D:/LOX1_EVAR_MR/results/tables/mediator_finngen_outcomes.csv", "utf-8");
const step1Csv = fs.readFileSync("D:/LOX1_EVAR_MR/results/tables/step1_final_summary.csv", "utf-8");
const exposureCsv = fs.readFileSync("D:/LOX1_EVAR_MR/data/exposure/eqtlgen_OLR1_blood_eqtl.csv", "utf-8");

function parseCsv(text) {
  const lines = text.trim().split("\n");
  const headers = lines[0].split(",");
  return lines.slice(1).map(line => {
    const vals = [];
    let cur = "", inQ = false;
    for (let i = 0; i < line.length; i++) {
      if (line[i] === '"') { inQ = !inQ; continue; }
      if (line[i] === "," && !inQ) { vals.push(cur); cur = ""; continue; }
      cur += line[i];
    }
    vals.push(cur);
    return vals;
  });
}

// Parse single SNP data
const singleSnpData = parseCsv(singleSnpCsv);
const looData = parseCsv(looCsv);
const medData = parseCsv(medOutcomesCsv);
const step1Data = parseCsv(step1Csv);
const exposureData = parseCsv(exposureCsv);

// Format p-value
function fmtP(p) {
  const v = parseFloat(p);
  if (v < 0.001) return v.toExponential(2);
  return v.toFixed(3);
}

// Format beta
function fmtB(b) {
  const v = parseFloat(b);
  return v.toFixed(4);
}

const doc = new Document({
  styles: {
    default: { document: { run: { font: "Times New Roman", size: 22 } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 26, bold: true, font: "Times New Roman" },
        paragraph: { spacing: { before: 240, after: 120 }, outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 24, bold: true, font: "Times New Roman" },
        paragraph: { spacing: { before: 200, after: 100 }, outlineLevel: 1 } },
    ],
  },
  sections: [{
    properties: { page: { size: { width: A4.width, height: A4.height }, margin: MARGINS } },
    footers: {
      default: new Footer({
        children: [new Paragraph({
          alignment: AlignmentType.CENTER,
          children: [
            new TextRun({ text: "Supplementary Materials \u2014 Page ", size: 18, font: "Times New Roman" }),
            new TextRun({ children: [PageNumber.CURRENT], size: 18, font: "Times New Roman" }),
          ],
        })],
      }),
    },
    children: [
      // Title
      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { before: 400, after: 120 },
        children: [new TextRun({ text: "Supplementary Materials", bold: true, size: 30, font: "Times New Roman" })],
      }),
      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { after: 240 },
        children: [new TextRun({ text: "Genetically Predicted OLR1 Expression and Abdominal Aortic Aneurysm Risk: A Mendelian Randomization Study with Colocalization and Mediation Analysis", italics: true, size: 22, font: "Times New Roman" })],
      }),
      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { after: 360 },
        children: [new TextRun({ text: "Wei Wei, Peter J. Holt, Yi-Kuan Chen", size: 20, font: "Times New Roman" })],
      }),

      // === Supplementary Table 1: Instrument Characteristics ===
      tableTitle("Supplementary Table 1. Characteristics of 31 OLR1 cis-eQTL SNPs from eQTLGen consortium."),
      body("SNPs are ordered by eQTL effect size (beta). Effect allele (EA) and other allele (OA) refer to the allele associated with OLR1 expression. Negative beta indicates decreased OLR1 expression. F-statistics were all >10 (not shown). 24 of 31 SNPs were successfully harmonized with FinnGen AAA GWAS and used in the primary analysis."),
      makeTable(
        ["SNP", "EA", "OA", "\u03b2 (eQTL)", "p (eQTL)"],
        exposureData.map(row => [row[0], row[3], row[4], fmtB(row[5]), fmtP(row[7])]),
        [2400, 1200, 1200, 2000, 2000]
      ),
      spacer(),

      // === Supplementary Table 2: Single-SNP Wald Ratios ===
      tableTitle("Supplementary Table 2. Single-SNP Wald ratio estimates for the eQTLGen primary analysis (24 SNPs)."),
      body("Each row represents the causal effect estimate (Wald ratio) when using only that single SNP as an instrument. The last three rows show the combined estimates from IVW, MR-Egger, and weighted median methods. \u03b2 = log odds ratio per standard deviation increase in OLR1 expression."),
      makeTable(
        ["SNP", "\u03b2 (SE)", "p-value"],
        singleSnpData.map(row => [row[5], `${fmtB(row[6])} (${parseFloat(row[7]).toFixed(5)})`, fmtP(row[8])]),
        [3000, 3500, 2500]
      ),
      spacer(),

      // === Supplementary Table 3: Leave-One-Out ===
      tableTitle("Supplementary Table 3. Leave-one-out sensitivity analysis (eQTLGen primary analysis)."),
      body("Each row shows the IVW estimate when the indicated SNP is removed from the analysis. The last row (All) shows the estimate using all 24 SNPs. Directional consistency across all iterations confirms that no single SNP drives the association."),
      makeTable(
        ["SNP excluded", "\u03b2 (SE)", "p-value"],
        looData.map(row => {
          const snpName = row[5].replace("All - Inverse variance weighted", "All (IVW)").replace("All - MR Egger", "All (MR Egger)");
          return [snpName, `${fmtB(row[6])} (${parseFloat(row[7]).toFixed(5)})`, fmtP(row[8])];
        }),
        [3000, 3500, 2500]
      ),
      spacer(),

      // === Supplementary Table 4: Heterogeneity and Pleiotropy ===
      tableTitle("Supplementary Table 4. Heterogeneity and pleiotropy tests for the eQTLGen primary analysis."),
      body("Cochran's Q test assesses heterogeneity across instrumental variable estimates. A non-significant p-value (p>0.05) indicates no evidence of heterogeneity. The MR-Egger intercept test assesses directional horizontal pleiotropy; a non-significant intercept (p>0.05) indicates no evidence of pleiotropy."),
      makeTable(
        ["Test", "Method", "Statistic", "df", "p-value", "Conclusion"],
        [
          ["Heterogeneity", "MR-Egger", "Q=2.39", "22", ">0.99", "No heterogeneity"],
          ["Heterogeneity", "IVW", "Q=2.91", "23", ">0.99", "No heterogeneity"],
          ["Pleiotropy", "MR-Egger intercept", "0.026 (SE=0.037)", "\u2014", "0.48", "No pleiotropy"],
        ],
        [1800, 2000, 1800, 800, 1200, 1400]
      ),
      spacer(),

      // === Supplementary Table 5: Mediator Outcomes ===
      tableTitle("Supplementary Table 5. Mediator instrument SNPs and their FinnGen AAA outcome data."),
      body("Effect estimates from FinnGen R12 AAA GWAS (I9_ABAORTANEUR, 4,439 cases / 463,106 controls). beta_aaa = log odds ratio for AAA per effect allele. CRP instruments are located on chr1 (CRP gene locus). IL6R instruments are located on chr1 (IL6R gene locus). rs2228145 (Asp358Ala) is the lead functional variant increasing sIL6R by ~34% per allele."),
      makeTable(
        ["Mediator", "SNP", "Chr", "Position", "EA/OA", "\u03b2_AAA (SE)", "p_AAA"],
        medData.map(row => [
          row[0], row[1], row[2], row[3],
          `${row[4]}/${row[5]}`,
          `${parseFloat(row[6]).toFixed(4)} (${parseFloat(row[7]).toFixed(4)})`,
          fmtP(row[8]),
        ]),
        [1200, 1600, 600, 1600, 1000, 2200, 1200]
      ),
      spacer(),

      // === Supplementary Table 6: Step 1 Summary ===
      tableTitle("Supplementary Table 6. Step 1 mediation MR summary: genetic overlap between OLR1 instruments and CRP/IL6R."),
      body("Step 1 tested whether OLR1 cis-eQTL SNPs (n=31) have significant trans-effects on CRP or IL6R at the mRNA level (eQTLGen trans-eQTL, 59,786 significant associations) or protein level (GWAS Catalog, 212 top CRP associations). The consistently null results across all approaches indicate that OLR1 and IL6R/CRP operate through independent parallel pathways rather than a linear cascade."),
      makeTable(
        ["Step", "Pathway", "Method", "Result", "p-value", "Conclusion"],
        step1Data.map(row => [row[0], row[1], row[2], row[3], row[4], row[5]]),
        [800, 2000, 1800, 1800, 1200, 1400]
      ),
      spacer(),

      // === Supplementary Table 7: Colocalization Summary ===
      tableTitle("Supplementary Table 7. Colocalization posterior probabilities across all three eQTL sources."),
      body("Coloc.abf Bayesian analysis testing five hypotheses: H0 (no association with either trait), H1 (association with eQTL only), H2 (association with AAA only), H3 (both associated, independent variants), H4 (both associated, shared causal variant). PP.H4 > 0.7 is considered evidence of colocalization. All three analyses show PP.H4 < 0.02, rejecting the shared causal variant hypothesis."),
      makeTable(
        ["eQTL Source", "nSNP", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "Conclusion"],
        [
          ["eQTLGen Blood", "31", "6\u00d710\u207b\u2076", "0.988", "3.2\u00d710\u207b\u2079", "0.0005", "0.011", "Reject H4"],
          ["Artery_Tibial", "49", "0.002", "0.979", "2.0\u00d710\u207b\u2076", "0.0008", "0.017", "Reject H4"],
          ["Artery_Coronary", "9", "0.235", "0.754", "3.3\u00d710\u207b\u2075", "0.0001", "0.011", "Reject H4"],
        ],
        [1600, 700, 1000, 800, 1200, 1000, 800, 926]
      ),
      spacer(),

      // === Supplementary Figure Legends ===
      heading("Supplementary Figure Legends", 1),

      tableTitle("Supplementary Figure 1. Scatter plot of SNP effects on OLR1 expression vs. AAA risk (eQTLGen)."),
      body("Each point represents a single SNP. The x-axis shows the effect of the SNP on OLR1 expression (beta from eQTLGen whole blood eQTL). The y-axis shows the effect of the SNP on AAA risk (log odds ratio from FinnGen R12). The slope of the line represents the IVW MR estimate. Positive slope indicates that higher OLR1 expression is associated with increased AAA risk. File: scatter_eQTLGen_Blood.png"),

      tableTitle("Supplementary Figure 2. Forest plot of single-SNP and combined MR estimates (eQTLGen)."),
      body("Each row shows the Wald ratio estimate (black dot) and 95% confidence interval (horizontal line) for an individual SNP. The bottom rows show combined estimates from IVW (diamond), MR-Egger, and weighted median methods. All individual SNP estimates are directionally consistent (positive). File: forest_eQTLGen_Blood.png"),

      tableTitle("Supplementary Figure 3. Funnel plot for asymmetry assessment (eQTLGen)."),
      body("Each point represents a single SNP's Wald ratio estimate plotted against its standard error. Symmetry of the funnel indicates no small-study bias or directional pleiotropy. The MR-Egger regression line (if shown) tests for funnel plot asymmetry. The symmetric distribution supports the absence of directional pleiotropy. File: funnel_eQTLGen_Blood.png"),

      tableTitle("Supplementary Figure 4. Leave-one-out forest plot (eQTLGen)."),
      body("Each row shows the IVW estimate (black dot) and 95% confidence interval when the indicated SNP is excluded from the analysis. The bottom row (All) shows the estimate using all 24 SNPs. All estimates remain positive and directionally consistent, confirming that no single SNP drives the association. File: loo_eQTLGen_Blood.png"),

      tableTitle("Supplementary Figure 5. Colocalization sensitivity plots."),
      body("Three panels showing colocalization results for (A) eQTLGen Blood, (B) GTEx Artery_Tibial, and (C) GTEx Artery_Coronary. Each panel displays: (top) the regional association plot for OLR1 eQTL signal, (middle) the regional association plot for AAA GWAS signal, (bottom) the posterior probabilities for H0-H4 at each SNP. In all three panels, PP.H4 (red) remains near zero across the locus, while PP.H1 (blue) dominates, indicating no shared causal variant. Files: coloc_sensitivity_*.png"),

      tableTitle("Supplementary Figure 6. Mediation MR forest plot (CRP and IL6R \u2192 AAA)."),
      body("Forest plot showing the causal effect estimates for CRP \u2192 AAA (top panel, 4 SNPs) and IL6R \u2192 AAA (bottom panel, 4 SNPs). Each row shows an individual SNP's Wald ratio, and the bottom rows show combined IVW, MR-Egger, and weighted median estimates. CRP shows no significant effect (all confidence intervals crossing null). IL6R shows a highly significant protective effect (all estimates well below null, OR\u22480.72). File: mediation_forest_plot.svg"),

      // === Main Figure Legends ===
      heading("Main Figure Legends", 1),

      tableTitle("Figure 1. Study design and analytical framework."),
      body("Schematic overview of the three-component Mendelian randomization study: (1) Total effect MR using OLR1 cis-eQTLs from eQTLGen (n=24 SNPs, whole blood) and GTEx v8 (Artery_Tibial, Artery_Coronary) as exposure instruments, with FinnGen R12 AAA GWAS (4,439 cases / 463,106 controls) as outcome. (2) Colocalization analysis (coloc.abf) testing whether eQTL and AAA signals share a causal variant at the OLR1 locus (\u00b1100 kb). (3) Mediation MR with two steps: Step 1 (OLR1\u2192mediator, trans-eQTL + GWAS Catalog) and Step 2 (mediator\u2192AAA, IVW MR). Key results are annotated at each step."),

      tableTitle("Figure 2. Integrated pathway model: parallel inflammatory pathways converging on AAA."),
      body("Schematic representation of the findings. Left panel: LOX-1/OLR1 pathway (chr12) \u2014 genetically upregulated OLR1 expression increases AAA risk (IVW OR=1.0023, p=0.022) through an indirect mechanism (coloc PP.H4 < 2%). Right panel: IL6R trans-signaling pathway (chr1) \u2014 increased sIL6R decreases AAA risk (IVW OR=0.722, p=2.96\u00d710\u207b\u00b2\u2076, 28% risk reduction). Center: CRP shows no causal effect on AAA (OR=0.968, p=0.58). Step 1 genetic overlap analysis (0/31 SNPs) confirms that OLR1 and IL6R operate through independent parallel pathways. The model supports dual-target intervention strategies (LOX-1 inhibition + IL6R modulation)."),
    ],
  }],
});

(async () => {
  const buffer = await Packer.toBuffer(doc);
  const outPath = "D:/workbuddy工作空间/LOX-1_EVAR_MR_Supplementary.docx";
  fs.writeFileSync(outPath, buffer);
  console.log(`Supplementary materials written: ${outPath} (${(buffer.length / 1024).toFixed(1)} KB)`);
})();

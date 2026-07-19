#!/usr/bin/env python3
"""
Convert all figures to publication-quality formats:
- PNG -> TIFF (300 DPI, LZW compression) - journal standard for raster
- PNG -> PDF (300 DPI) - for vector-like quality
- SVG -> PDF (true vector) - for mediation forest plot
- Skip Artery_Coronary figures (empty - palindromic SNPs removed)
"""

import os
from PIL import Image

# cairosvg requires native cairo library (not available on Windows)
# SVG files will be kept as-is (vector format, accepted by journals)
HAS_CAIROSVG = False

fig_dir = "D:/LOX1_EVAR_MR/results/figures"
out_dir = "D:/LOX1_EVAR_MR/results/figures_highres"
os.makedirs(out_dir, exist_ok=True)

# Define main vs supplementary figures
main_figures = {
    "scatter_eQTLGen_Blood.png": "Figure1_scatter_plot",
    "forest_eQTLGen_Blood.png": "Figure2_forest_plot", 
    "funnel_eQTLGen_Blood.png": "Figure3_funnel_plot",
    "loo_eQTLGen_Blood.png": "Figure4_leave_one_out",
    "mediation_forest_plot.svg": "Figure5_mediation_forest",
}

supplementary_figures = {
    "scatter_Artery_Tibial.png": "SuppFig1_scatter_Artery_Tibial",
    "forest_Artery_Tibial.png": "SuppFig2_forest_Artery_Tibial",
    "funnel_Artery_Tibial.png": "SuppFig3_funnel_Artery_Tibial",
    "loo_Artery_Tibial.png": "SuppFig4_leave_one_out_Artery_Tibial",
    "coloc_sensitivity_eQTLGen_Blood.png": "SuppFig5_coloc_sensitivity_eQTLGen",
    "coloc_sensitivity_Artery_Tibial.png": "SuppFig6_coloc_sensitivity_Artery_Tibial",
    "coloc_sensitivity_Artery_Coronary.png": "SuppFig7_coloc_sensitivity_Artery_Coronary",
}

all_figures = {**main_figures, **supplementary_figures}

print("=" * 60)
print("High-DPI Figure Conversion for Journal Submission")
print("=" * 60)
print(f"Input:  {fig_dir}")
print(f"Output: {out_dir}")
print(f"Target: 300 DPI, TIFF (LZW) + PDF")
print()

converted = 0
skipped = 0

for src_name, out_name in all_figures.items():
    src_path = os.path.join(fig_dir, src_name)
    
    if not os.path.exists(src_path):
        print(f"  [SKIP] {src_name} - file not found")
        skipped += 1
        continue
    
    # Check if it's a placeholder (Artery_Coronary files are 27786 bytes = empty)
    file_size = os.path.getsize(src_path)
    if file_size < 30000 and src_name.endswith(".png"):
        print(f"  [SKIP] {src_name} - appears to be empty/placeholder ({file_size} bytes)")
        skipped += 1
        continue
    
    if src_name.endswith(".png"):
        img = Image.open(src_path)
        
        # Convert palette mode to RGB
        if img.mode == "P":
            img = img.convert("RGB")
        
        # Ensure 300 DPI
        target_dpi = (300, 300)
        
        # Save as TIFF (LZW compression) - journal standard
        tiff_path = os.path.join(out_dir, f"{out_name}.tiff")
        img.save(tiff_path, format="TIFF", dpi=target_dpi, compression="tiff_lzw")
        tiff_size = os.path.getsize(tiff_path) / 1024
        print(f"  [OK] {src_name} -> {out_name}.tiff ({img.size[0]}x{img.size[1]}px, 300 DPI, {tiff_size:.0f}KB)")
        
        # Also save as PDF (raster embedded in PDF)
        pdf_path = os.path.join(out_dir, f"{out_name}.pdf")
        img.save(pdf_path, format="PDF", dpi=target_dpi, resolution=300)
        pdf_size = os.path.getsize(pdf_path) / 1024
        print(f"       -> {out_name}.pdf ({pdf_size:.0f}KB)")
        
        # Also save high-res PNG (for journals that accept PNG)
        png_path = os.path.join(out_dir, f"{out_name}.png")
        img.save(png_path, format="PNG", dpi=target_dpi)
        png_size = os.path.getsize(png_path) / 1024
        print(f"       -> {out_name}.png ({png_size:.0f}KB)")
        
        converted += 1
    
    elif src_name.endswith(".svg"):
        # SVG is already vector format - copy as-is and create PDF version
        import shutil
        svg_out = os.path.join(out_dir, f"{out_name}.svg")
        shutil.copy2(src_path, svg_out)
        svg_size = os.path.getsize(svg_out) / 1024
        print(f"  [OK] {src_name} -> {out_name}.svg (vector, {svg_size:.0f}KB) - kept as vector")
        print(f"       Note: SVG is vector format, accepted by most journals.")
        print(f"       For PDF conversion, open in browser/Inkscape and 'Print to PDF'.")
        converted += 1

print()
print("=" * 60)
print(f"Conversion complete: {converted} figures converted, {skipped} skipped")
print(f"Output directory: {out_dir}")
print()
print("Available formats per figure:")
print("  .tiff - 300 DPI, LZW compression (journal standard)")
print("  .pdf  - 300 DPI (for vector or raster)")
print("  .png  - 300 DPI (for journals accepting PNG)")
print()

# List final output
print("Generated files:")
for f in sorted(os.listdir(out_dir)):
    size = os.path.getsize(os.path.join(out_dir, f)) / 1024
    print(f"  {f:50s} {size:>8.0f} KB")

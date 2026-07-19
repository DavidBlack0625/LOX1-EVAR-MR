#!/usr/bin/env python3
"""
Render the mediation forest plot as high-res TIFF/PNG (300 DPI)
using Pillow, since cairosvg requires native cairo library.
"""
from PIL import Image, ImageDraw, ImageFont
import math

# Canvas: 300 DPI, ~8x6 inches
W, H = 2400, 1600
SCALE = W / 900  # Scale factor from SVG coordinates

img = Image.new("RGB", (W, H), "#FAFAFA")
draw = ImageDraw.Draw(img)

# Try to load a TrueType font
try:
    font_title = ImageFont.truetype("C:/Windows/Fonts/arialbd.ttf", int(18 * SCALE))
    font_sub = ImageFont.truetype("C:/Windows/Fonts/arial.ttf", int(12 * SCALE))
    font_bold = ImageFont.truetype("C:/Windows/Fonts/arialbd.ttf", int(11 * SCALE))
    font_normal = ImageFont.truetype("C:/Windows/Fonts/arial.ttf", int(10 * SCALE))
    font_small = ImageFont.truetype("C:/Windows/Fonts/arial.ttf", int(9 * SCALE))
except:
    font_title = ImageFont.load_default()
    font_sub = ImageFont.load_default()
    font_bold = ImageFont.load_default()
    font_normal = ImageFont.load_default()
    font_small = ImageFont.load_default()

def sx(x): return int(x * SCALE)
def sy(y): return int(y * SCALE)
def ss(s): return int(s * SCALE)

# Title
draw.text((sx(450), sy(15)), "Mediation MR: Forest Plot", fill="#333",
          font=font_title, anchor="mm")
draw.text((sx(450), sy(42)), "OLR1 -> Mediator -> AAA (FinnGen R12, I9_ABAORTANEUR)",
          fill="#666", font=font_sub, anchor="mm")

# X-axis
x_axis_y = sy(505)
x_axis_left = sx(250)
x_axis_right = sx(700)
draw.line([(x_axis_left, x_axis_y), (x_axis_right, x_axis_y)], fill="#333", width=ss(2))

# X-axis ticks and labels
ticks = [0.67, 0.74, 0.82, 0.90, 1.00, 1.11, 1.22, 1.35]
for tick_val in ticks:
    # Map OR to x position: log scale, OR=1.0 at x=531.25
    log_val = math.log(tick_val)
    log_1 = math.log(1.0)
    # In SVG: x = 531.25 + log_val * (475.0 - 531.25) / log(0.90)
    # Simplify: scale = (531.25 - 475.0) / |log(0.90)| = 56.25 / 0.10536 = 533.8
    x_pos = sx(531.25 + log_val * 533.8)
    draw.line([(x_pos, x_axis_y), (x_pos, x_axis_y + ss(5))], fill="#666", width=ss(1))
    draw.text((x_pos, sy(520)), f"{tick_val:.2f}", fill="#666", font=font_small, anchor="mm")

# OR=1.0 reference line
ref_x = sx(531.25)
draw.line([(ref_x, sy(65)), (ref_x, sy(505))], fill="#999", width=ss(1))
# Dashed line effect
for y in range(sy(65), sy(505), ss(7)):
    draw.line([(ref_x, y), (ref_x, min(y + ss(4), sy(505)))], fill="#999", width=ss(1))

draw.text((ref_x, sy(58)), "OR=1.0", fill="#999", font=font_small, anchor="mm")
draw.text((sx(475), sy(540)), "Odds Ratio (95% CI)", fill="#333", font=font_normal, anchor="mm")

# Column headers
draw.text((sx(15), sy(78)), "Analysis", fill="#333", font=font_bold)
draw.text((sx(710), sy(78)), "OR (95% CI)", fill="#333", font=font_bold)
draw.text((sx(815), sy(78)), "p-value", fill="#333", font=font_bold)

def draw_row(y_svg, label, or_val, ci_lo, ci_hi, p_text, color, is_bold=False, snp_text=""):
    """Draw a forest plot row."""
    y = sy(y_svg)
    # Label
    draw.text((sx(15), y), label, fill="#333", font=font_bold if is_bold else font_normal)
    if snp_text:
        draw.text((sx(15), y + ss(14)), snp_text, fill="#666", font=font_small)
    
    # Calculate x positions on log scale
    def or_to_x(or_v):
        return sx(531.25 + math.log(or_v) * 533.8)
    
    x_lo = or_to_x(ci_lo)
    x_hi = or_to_x(ci_hi)
    x_pt = or_to_x(or_val)
    
    # CI line
    draw.line([(x_lo, y), (x_hi, y)], fill=color, width=ss(2))
    # End caps
    draw.line([(x_lo, y - ss(4)), (x_lo, y + ss(4))], fill=color, width=ss(2))
    draw.line([(x_hi, y - ss(4)), (x_hi, y + ss(4))], fill=color, width=ss(2))
    # Point estimate (square)
    sq_size = ss(5) if not is_bold else ss(5)
    draw.rectangle([(x_pt - sq_size, y - sq_size), (x_pt + sq_size, y + sq_size)],
                   fill=color)
    
    # Text columns
    ci_text = f"{or_val:.3f} ({ci_lo:.3f}-{ci_hi:.3f})"
    draw.text((sx(710), y), ci_text, fill="#333", font=font_normal)
    draw.text((sx(815), y), p_text, fill=color, font=font_bold)

def draw_section_header(y_svg, title):
    """Draw a section header with background."""
    y = sy(y_svg)
    draw.rectangle([(sx(5), y), (sx(895), y + ss(22))], fill="#E8EAF6")
    draw.text((sx(15), y + ss(6)), title, fill="#333", font=font_bold)

# Total Effect section
draw_section_header(95, "Total Effect")
draw_row(132, "IVW (eQTLGen, 24 SNPs)", 1.0023, 1.0003, 1.0043, "0.022 *", "#2196F3")

# CRP section
draw_section_header(153, "CRP (Step 2)")
draw_row(190, "Inverse variance weighted (4 SNPs)", 0.968, 0.864, 1.085, "0.580", "#9E9E9E")
draw_row(220, "MR Egger (4 SNPs)", 1.013, 0.042, 24.357, "0.994", "#9E9E9E")
draw_row(250, "Weighted median (4 SNPs)", 0.927, 0.828, 1.039, "0.192", "#9E9E9E")

# IL6R section
draw_section_header(271, "IL6R (Step 2)")
draw_row(308, "Inverse variance weighted (4 SNPs)", 0.722, 0.679, 0.767, "2.96e-26 ***", "#F44336")
draw_row(338, "MR Egger (4 SNPs)", 0.724, 0.682, 0.769, "5.25e-26 ***", "#F44336")
draw_row(368, "Weighted median (4 SNPs)", 0.714, 0.677, 0.752, "1.10e-36 ***", "#F44336")

# Save in multiple formats
out_dir = "D:/LOX1_EVAR_MR/results/figures_highres"
dpi = (300, 300)

# TIFF (LZW)
tiff_path = f"{out_dir}/Figure5_mediation_forest.tiff"
img.save(tiff_path, format="TIFF", dpi=dpi, compression="tiff_lzw")
print(f"Saved: {tiff_path} ({W}x{H}px, 300 DPI)")

# PNG
png_path = f"{out_dir}/Figure5_mediation_forest.png"
img.save(png_path, format="PNG", dpi=dpi)
print(f"Saved: {png_path}")

# PDF
pdf_path = f"{out_dir}/Figure5_mediation_forest.pdf"
img.save(pdf_path, format="PDF", dpi=dpi, resolution=300)
print(f"Saved: {pdf_path}")

print("\nMediation forest plot rendered successfully at 300 DPI!")

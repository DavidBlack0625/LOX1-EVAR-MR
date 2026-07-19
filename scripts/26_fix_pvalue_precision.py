"""
Fix p-value precision in mediation MR results.
The original norm_sf used 2*(1-nd.cdf(abs(x))) which loses precision for z>8.
Use math.erfc for high precision: p = erfc(|z| / sqrt(2))
"""
import csv, math, os

RESULTS_DIR = "D:/LOX1_EVAR_MR/results/tables"

def precise_p(z):
    """High-precision two-tailed p-value using erfc."""
    return math.erfc(abs(z) / math.sqrt(2))

# Read current results
results = []
with open(f"{RESULTS_DIR}/mediation_mr_results.csv", 'r') as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    for row in reader:
        results.append(row)

# Recalculate p-values
print("=== Fixing p-value precision ===")
for row in results:
    if row['pval'] == '0.0' or float(row['pval']) == 0.0:
        b = float(row['b'])
        se = float(row['se'])
        if se > 0:
            z = b / se
            new_p = precise_p(z)
            row['pval'] = new_p
            print(f"  {row['mediator']} {row['method']}: z={z:.4f}, p={new_p:.4e}")
        else:
            print(f"  {row['mediator']} {row['method']}: SE=0, cannot calculate p")
    else:
        old_p = float(row['pval'])
        if old_p > 0 and old_p < 1:
            # Also verify non-zero p-values are precise
            b = float(row['b'])
            se = float(row['se'])
            if se > 0:
                z = b / se
                new_p = precise_p(z)
                if abs(new_p - old_p) / old_p > 0.01:  # >1% difference
                    print(f"  {row['mediator']} {row['method']}: old_p={old_p:.4e} -> new_p={new_p:.4e} (corrected)")
                    row['pval'] = new_p
                else:
                    print(f"  {row['mediator']} {row['method']}: p={old_p:.4e} (OK)")

# Write fixed results
with open(f"{RESULTS_DIR}/mediation_mr_results.csv", 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for row in results:
        writer.writerow(row)
print(f"\nSaved fixed: {RESULTS_DIR}/mediation_mr_results.csv")

# Also fix summary table
summary_rows = []
with open(f"{RESULTS_DIR}/mediation_mr_summary.csv", 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        summary_rows.append(row)

print("\n=== Fixing summary table ===")
for i, row in enumerate(summary_rows):
    if i == 0:  # header
        continue
    pval_str = row[10]  # p-value column
    try:
        pval = float(pval_str)
        if pval == 0.0:
            # Need to recalculate from beta and SE
            beta = float(row[5])
            se = float(row[6])
            if se > 0:
                z = beta / se
                new_p = precise_p(z)
                row[10] = f"{new_p:.4e}"
                sig = '***' if new_p < 0.001 else ('**' if new_p < 0.01 else ('*' if new_p < 0.05 else 'ns'))
                row[11] = sig
                print(f"  {row[0]} {row[1]} {row[3]}: z={z:.4f}, p={new_p:.4e} {sig}")
    except (ValueError, IndexError):
        pass

with open(f"{RESULTS_DIR}/mediation_mr_summary.csv", 'w', newline='') as f:
    writer = csv.writer(f)
    for row in summary_rows:
        writer.writerow(row)
print(f"Saved fixed: {RESULTS_DIR}/mediation_mr_summary.csv")

# Print final summary
print("\n" + "=" * 70)
print("FINAL MEDIATION MR RESULTS (precision-fixed)")
print("=" * 70)
print(f"{'Step':<12} {'Path':<25} {'Method':<25} {'OR':<10} {'95% CI':<20} {'p-value':<12} {'Sig':<5}")
print("-" * 110)
for row in results:
    or_val = float(row['OR'])
    or_lo = float(row['OR_lower'])
    or_hi = float(row['OR_upper'])
    p = float(row['pval'])
    sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else 'ns'))
    print(f"{row['step']:<12} {row['mediator']:<25} {row['method']:<25} {or_val:.4f}    ({or_lo:.4f}-{or_hi:.4f})  {p:.4e}  {sig}")

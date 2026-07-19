import gzip, csv

exp_snps = set()
with open('D:/LOX1_EVAR_MR/data/exposure/eqtlgen_OLR1_blood_eqtl.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        exp_snps.add(row['SNP'])

print(f'eQTLGen OLR1 SNPs to find: {len(exp_snps)}')

matches = []
found = set()
finngen_gz = 'D:/workbuddy工作空间/LOX-1_EVAR_MR/data/outcome/finngen_R12_I9_ABAORTANEUR.gz'
with gzip.open(finngen_gz, 'rt') as f:
    header = f.readline().strip().split('\t')
    for line in f:
        parts = line.strip().split('\t')
        for rsid in parts[4].split(';'):
            rsid = rsid.strip()
            if rsid in exp_snps:
                row = dict(zip(header, parts))
                row['matched_snp'] = rsid
                matches.append(row)
                found.add(rsid)

print(f'Found {len(matches)} matches for {len(found)} unique SNPs')

if matches:
    with open('D:/LOX1_EVAR_MR/data/outcome/finngen_aaa_eqtlgen_matched.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['matched_snp'] + header)
        writer.writeheader()
        writer.writerows(matches)
    print('Saved to finngen_aaa_eqtlgen_matched.csv')
    for m in matches[:5]:
        print(f"  {m['matched_snp']}  beta={m['beta']}  p={m['pval']}")

import csv
from collections import defaultdict

entropy_path = '../Analysis_output/VEEV_entropy_genome_DPI3_DPI5.csv'
mut_path = '../Analysis_output/VEEV_LoFreq_mutations.csv'

bases = ['A', 'C', 'G', 'T']
acc = defaultdict(lambda: {'sumA': 0.0, 'sumC': 0.0, 'sumG': 0.0, 'sumT': 0.0, 'n': 0})


def parse_af(v):
    s = (v or '').strip()
    if not s:
        return []
    if ',' in s:
        out = []
        for t in s.split(','):
            t = t.strip()
            if not t:
                continue
            try:
                out.append(float(t))
            except Exception:
                pass
        return out
    try:
        return [float(s)]
    except Exception:
        return []


def norm_base(b):
    b = (b or '').strip().upper()
    return b if b in bases else None


with open(mut_path, newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            dpi = int(float(row['DPI']))
            pos = int(float(row['Position']))
        except Exception:
            continue
        if dpi not in (3, 5):
            continue

        ref = norm_base(row.get('Reference', ''))
        alts = [norm_base(x) for x in (row.get('Alternate', '') or '').split(',') if x.strip()]
        afs = parse_af(row.get('Allele_Frequency', ''))

        if not alts:
            continue
        if len(afs) < len(alts):
            afs = afs + [0.0] * (len(alts) - len(afs))
        elif len(afs) > len(alts):
            afs = afs[:len(alts)]
        if not afs:
            continue

        if max(afs) > 1.0:
            afs = [x / 100.0 for x in afs]
        afs = [min(1.0, max(0.0, x)) for x in afs]

        p_alt = min(1.0, sum(afs))
        p_ref = max(0.0, 1.0 - p_alt)

        prof = {b: 0.0 for b in bases}
        if ref:
            prof[ref] += p_ref
        for b, p in zip(alts, afs):
            if b:
                prof[b] += p

        s = sum(prof.values())
        if s <= 0:
            continue
        for b in bases:
            prof[b] /= s

        key = (dpi, pos)
        acc[key]['sumA'] += prof['A'] * 100.0
        acc[key]['sumC'] += prof['C'] * 100.0
        acc[key]['sumG'] += prof['G'] * 100.0
        acc[key]['sumT'] += prof['T'] * 100.0
        acc[key]['n'] += 1

rows = []
with open(entropy_path, newline='') as f:
    reader = csv.DictReader(f)
    fieldnames = list(reader.fieldnames)
    for c in ['BasePct_A', 'BasePct_C', 'BasePct_G', 'BasePct_T']:
        if c not in fieldnames:
            fieldnames.append(c)

    for row in reader:
        try:
            key = (int(float(row['DPI'])), int(float(row['Position'])))
        except Exception:
            key = None

        if key in acc and acc[key]['n'] > 0:
            n = acc[key]['n']
            row['BasePct_A'] = f"{acc[key]['sumA'] / n:.6f}"
            row['BasePct_C'] = f"{acc[key]['sumC'] / n:.6f}"
            row['BasePct_G'] = f"{acc[key]['sumG'] / n:.6f}"
            row['BasePct_T'] = f"{acc[key]['sumT'] / n:.6f}"
        else:
            row['BasePct_A'] = '0.000000'
            row['BasePct_C'] = '0.000000'
            row['BasePct_G'] = '0.000000'
            row['BasePct_T'] = '0.000000'
        rows.append(row)

with open(entropy_path, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f'Updated {entropy_path} with base percentages; rows={len(rows)}')

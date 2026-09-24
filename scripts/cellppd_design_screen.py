#!/usr/bin/env python3
"""CellPPD (Gautam et al.) per-sequence design endpoint screen of the 18 designed
CPP candidates. The batch prediction endpoint of CellPPD accepted jobs but never
populated results (dead backend, documented honest negative); the live path is
pep_test.php, which returns the full single-mutant landscape (all 19 single
substitutions x L positions) with CellPPD SVM+motif scores. We fetch the
landscape for each candidate and summarize mutational sensitivity.

Resumable: raw HTML cached in --rawdir; parsed summary appended per candidate.
"""
import argparse, json, os, re, subprocess, sys, time

def fetch(url, data_fields, timeout=60):
    cmd = ['curl', '-s', '-m', str(timeout)]
    if data_fields:
        cmd += ['-X', 'POST', url]
        for f in data_fields:
            cmd += ['-F', f]
    else:
        cmd += [url]
    return subprocess.run(cmd, capture_output=True, text=True).stdout

def submit(seq):
    out = fetch('https://webs.iiitd.edu.in/raghava/cellppd/pep_test.php',
        [f'seq={seq}', 'method=2', 'eval=10', 'thval=0.0',
         'field[]=4', 'field[]=7', 'field[]=9', 'field[]=11', 'field[]=13'])
    m = re.findall(r'url=(submitfreq\.php\?ran=\d+)', out)
    return m[0].split('ran=')[1] if m else None

def poll(ran, tries=8, sleep_s=12):
    html = ''
    for _ in range(tries):
        time.sleep(sleep_s)
        html = fetch(f'https://webs.iiitd.edu.in/raghava/cellppd/submitfreq.php?ran={ran}', None, 40)
        if len(html) > 10000:
            break
    return html

def parse_landscape(html, seq):
    """Rows in the (malformed) markup look like:
    <tr><font color='red'>A</font>RWRKKKRKK</a></td><td align=center>1</td><td align=center>0.00</td><td align=center>-1.02</td></tr>
    => mutant_aa, mutant_peptide, position, motif_indicator, svm_score"""
    rows = []
    pat = re.compile(
        r"<tr>([A-Z]*)<font color='red'>([A-Z])</font>([A-Z]*)</a></td>"
        r"<td align=center>(\d+)</td>"
        r"<td align=center>([-\d.]+)</td>"
        r"<td align=center>([-\d.]+)</td></tr>")
    for m in pat.finditer(html):
        rows.append({'mutant': m.group(2),
                     'mutant_peptide': m.group(1) + m.group(2) + m.group(3),
                     'position': int(m.group(4)), 'motif_indicator': float(m.group(5)),
                     'svm_score': float(m.group(6))})
    L = len(seq)
    if not rows:
        return None
    ok = all(sum(1 for r in rows if r['position'] == p) == 19 for p in range(1, L+1))
    scores = [r['svm_score'] for r in rows]
    per_pos = {}
    for p in range(1, L+1):
        ps = [r['svm_score'] for r in rows if r['position'] == p]
        if ps:
            per_pos[p] = {'mean': sum(ps)/len(ps), 'min': min(ps), 'max': max(ps)}
    worst_pos = max(per_pos, key=lambda p: per_pos[p]['max'] - per_pos[p]['min']) if per_pos else None
    return {
        'n_mutants': len(rows), 'length': L, 'positions_complete': ok,
        'svm_mean': sum(scores)/len(scores), 'svm_min': min(scores), 'svm_max': max(scores),
        'frac_mutants_cpp_at_thr0': sum(1 for s in scores if s > 0.0)/len(scores),
        'most_sensitive_position': worst_pos,
        'per_position': {str(p): v for p, v in sorted(per_pos.items())},
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--rawdir', default='/tmp/cellppd_raw')
    ap.add_argument('--out', default='results/cpp_cellppd_design.json')
    ap.add_argument('--max-candidates', type=int, default=99)
    args = ap.parse_args()
    os.makedirs(args.rawdir, exist_ok=True)

    cands = json.load(open('results/cpp_novel_candidates.json'))['named_candidates']
    items = list(cands.items()) if isinstance(cands, dict) else [(c['name'], c['sequence']) for c in cands]

    summary = {}
    if os.path.exists(args.out):
        summary = json.load(open(args.out)).get('candidates', {})

    done = 0
    for name, seq in items:
        if name in summary and summary[name].get('n_mutants'):
            continue
        if done >= args.max_candidates:
            break
        raw = os.path.join(args.rawdir, f'{name}.html')
        if os.path.exists(raw) and os.path.getsize(raw) > 10000:
            html = open(raw).read()
            src = 'cache'
        else:
            ran = submit(seq)
            if not ran:
                print(f'{name}: submission failed'); continue
            html = poll(ran)
            open(raw, 'w').write(html)
            src = f'ran={ran}'
            time.sleep(2)
        parsed = parse_landscape(html, seq)
        if parsed:
            parsed['sequence'] = seq
            summary[name] = parsed
            print(f'{name}: OK ({src}) n={parsed["n_mutants"]} mean={parsed["svm_mean"]:.3f} complete={parsed["positions_complete"]}')
        else:
            print(f'{name}: PARSE FAIL ({src}) bytes={len(html)}')
        done += 1

    out = {
        'tool': 'CellPPD (SVM+motif, Gautam et al., webs.iiitd.edu.in/raghava/cellppd)',
        'endpoint': 'pep_test.php single-sequence design (batch endpoint multiple_test.php accepted jobs but never returned rows - documented dead backend, not used for claims)',
        'method': 'SVM+motif (method=2), E-value 10, threshold 0.0',
        'date': time.strftime('%Y-%m-%d'),
        'candidates': summary,
    }
    json.dump(out, open(args.out, 'w'), indent=1)
    print(f'Wrote {args.out} with {len(summary)} candidates')

if __name__ == '__main__':
    main()

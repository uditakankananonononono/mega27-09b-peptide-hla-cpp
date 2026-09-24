"""ToxinPred batch toxicity screen of the 18 designed CPPs (tool: ToxinPred, IIITD Raghava Lab).

Submits the 18 novel CPP candidates (results/cpp_novel_candidates.json) to the
ToxinPred batch endpoint (SVM (Swiss-Prot) + Motif model, threshold 0.0) and
parses the per-peptide SVM scores and toxicity calls into
results/cpp_toxinpred_screen.json.

The server flow is: POST multipart form to multiple_test.php -> meta-refresh to
multi_submitfreq_S_motif.php?ran=<id> -> poll until the results table renders.

Run: python3 scripts/toxinpred_screen.py   (requires network access)
"""
import json, re, subprocess, sys, time, urllib.request

BASE = 'https://webs.iiitd.edu.in/raghava/toxinpred/'


def submit(sequences):
    seq_block = '\n'.join(sequences)
    cmd = ['curl', '-s', '-c', '/tmp/tp_cookies.txt', '-b', '/tmp/tp_cookies.txt',
           '-m', '120', '-A', 'Mozilla/5.0', '-F', f'seq={seq_block}',
           '-F', 'method=2', '-F', 'eval=10', '-F', 'thval=0.0',
           '-F', 'field[]=4', '-F', 'field[]=7', '-F', 'field[]=9',
           '-F', 'field[]=11', '-F', 'field[]=13', BASE + 'multiple_test.php']
    html = subprocess.run(cmd, capture_output=True, text=True).stdout
    m = re.search(r"url=(multi_submitfreq_S_motif\.php\?ran=\d+)", html)
    if not m:
        raise RuntimeError('no results redirect in submission response')
    return BASE + m.group(1)


def poll(url, tries=10, delay=20):
    for _ in range(tries):
        html = urllib.request.urlopen(urllib.request.Request(
            url, headers={'User-Agent': 'Mozilla/5.0'}), timeout=90).read().decode('latin-1')
        if 'http-equiv=\'refresh\'' not in html and html.count('Non-Toxin') + html.count('>Toxin<') > 0:
            return html
        time.sleep(delay)
    raise RuntimeError('ToxinPred job did not finish in time')


def parse(html):
    rows = []
    for r in re.findall(r'<tr[^>]*>(.*?)</tr>', html, re.S):
        tds = [re.sub(r'<[^>]*>', '', t).strip() for t in re.findall(r'<td[^>]*>(.*?)</td>', r, re.S)]
        if len(tds) >= 4 and any(re.fullmatch(r'[A-Z]{5,}', t) for t in tds):
            rows.append(tds)
    return rows


def main():
    cands = json.load(open('results/cpp_novel_candidates.json'))['named_candidates']
    by_seq = {c['sequence']: c['name'] for c in cands}
    html = poll(submit(list(by_seq)))
    out = {'tool': 'ToxinPred (IIITD Raghava Lab), batch submission, SVM (Swiss-Prot) + Motif model, threshold 0.0',
           'tool_url': BASE + 'multi_submit.php', 'results': []}
    for _, seq, svm, pred, *rest in parse(html):
        out['results'].append({'name': by_seq.get(seq, '?'), 'sequence': seq,
                               'svm_score': float(svm), 'prediction': pred})
    json.dump(out, open('results/cpp_toxinpred_screen.json', 'w'), indent=1)
    print(f"{sum(1 for r in out['results'] if r['prediction'] == 'Non-Toxin')}/{len(out['results'])} Non-Toxin")


if __name__ == '__main__':
    sys.exit(main())

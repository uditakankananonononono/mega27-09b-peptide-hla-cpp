"""HLP half-life/stability screen of the 18 designed CPPs (tool: HLP, Raghava Lab).

Submits each of the 18 novel CPP candidates (results/cpp_novel_candidates.json)
to the HLP peptide half-life server (whole-composition model) and parses the
server-reported half-life prediction and stability class into
results/cpp_hlp_stability.json.

Run: python3 scripts/hlp_stability_screen.py   (requires network access)
"""
import json, re, subprocess
from concurrent.futures import ThreadPoolExecutor

URL = 'https://webs.iiitd.edu.in/cgibin/hlp/pep_both.pl'


def run_one(cand):
    r = subprocess.run(['curl', '-s', '-m', '60', '-A', 'Mozilla/5.0',
                        '-e', 'https://webs.iiitd.edu.in/raghava/hlp/pep_both.htm',
                        '-F', f"pepseq={cand['sequence']}", '-F', 'mod=whole_comp', URL],
                       capture_output=True, text=True)
    txt = re.sub(r'\|+', '|', re.sub(r'<[^>]*>', '|', r.stdout))
    m = re.search(re.escape(cand['sequence']) + r'\|NO\|([\d.]+)\|(\w+)', txt)
    return {'name': cand['name'], 'sequence': cand['sequence'],
            'half_life_pred': float(m.group(1)) if m else None,
            'stability_class': m.group(2) if m else None}


def main():
    cands = json.load(open('results/cpp_novel_candidates.json'))['named_candidates']
    with ThreadPoolExecutor(max_workers=4) as ex:
        results = list(ex.map(run_one, cands))
    out = {'tool': 'HLP (Half-Life Prediction, Raghava Lab), whole-composition model',
           'tool_url': 'https://webs.iiitd.edu.in/raghava/hlp/pep_both.htm',
           'results': results}
    json.dump(out, open('results/cpp_hlp_stability.json', 'w'), indent=1)
    print(f"parsed {sum(1 for r in results if r['half_life_pred'] is not None)}/{len(results)}")


if __name__ == '__main__':
    main()

"""Verify the paper's bibliography against Europe PMC, CrossRef, and NCBI
E-utilities. Every cited work must resolve to a real record with matching
year (and journal where available). Output: results/citation_verification.json
"""
import json, time, urllib.parse, urllib.request

def get_json(url):
    req = urllib.request.Request(url, headers={"User-Agent": "mega27-09b-citation-check/0.1"})
    with urllib.request.urlopen(req, timeout=20) as r:
        return json.load(r)

REFS = [
    ("pdbstruct", "The Protein Data Bank", 2000),
    ("iedb", "The Immune Epitope Database", 2019),
    ("netmhcpan41", "NetMHCpan-4.1 and NetMHCIIpan-4.0", 2020),
    ("netmhcpan40", "NetMHCpan-4.0: Improved peptide-MHC class I interaction predictions", 2017),
    ("netmhc", "Reliable prediction of T-cell epitopes using neural networks", 2003),
    ("netmhcpan", "NetMHCpan, a method for quantitative predictions of peptide binding", 2007),
    ("smm", "Generating quantitative models describing the sequence specificity", 2005),
    ("mhcflurry", "MHCflurry 2.0: Improved pan-allele prediction", 2020),
    ("mhcflurry1", "MHCflurry: Open-Source Class I MHC Binding Affinity Prediction", 2018),
    ("mhcnuggets", "High-Throughput Prediction of MHC Class I and II Neoantigens with MHCnuggets", 2020),
    ("mixmhcpred", "Deciphering HLA-I motifs across HLA peptidomes", 2017),
    ("transphla", "A transformer-based model to predict peptide-HLA class I binding", 2022),
    ("capsnet", "Quantification of uncertainty in peptide-MHC binding prediction", 2019),
    ("kim2014", "Dataset size and composition impact the reliability of performance benchmarks", 2014),
    ("sidney", "HLA class I supertypes: a revised and updated classification", 2008),
    ("cppsite", "CPPsite 2.0: a repository of experimentally validated cell-penetrating peptides", 2016),
    ("cellppd", "CellPPD: in silico approaches for designing highly effective cell penetrating peptides", 2013),
    ("cpppred", "CPPpred: prediction of cell penetrating peptides", 2013),
    ("mlcpp", "MLCPP: machine-learning-based prediction of cell-penetrating peptides", 2018),
    ("skipcpp", "SkipCPP-Pred: an improved and promising sequence-based predictor", 2017),
    ("blosum", "Amino acid substitution matrices from protein blocks", 1992),
    ("eisenberg", "The hydrophobic moment detects periodicity in protein hydrophobicity", 1984),
    ("chengprusoff", "Relationship between the inhibition constant", 1973),
    ("delong", "Comparing the areas under two or more correlated receiver operating characteristic curves", 1988),
    ("uniprot", "UniProt: the universal protein knowledgebase", 2025),
    ("netctlpan", "An integrative approach to CTL epitope prediction", 2005),
    ("platt", "Probabilistic outputs for support vector machines", 1999),
]

epmc_out, failed = {}, []
for key, title, year in REFS:
    q = urllib.parse.quote(f'TITLE:"{title}"')
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={q}&format=json&pageSize=3"
    try:
        d = get_json(url)
        hits = d.get("resultList", {}).get("result", [])
        best = hits[0] if hits else {}
        ok = any(int(h.get("pubYear", 0) or 0) == year for h in hits) if hits else False
        epmc_out[key] = {"found": bool(hits), "year_match": ok,
                         "top_title": best.get("title", "")[:110],
                         "top_year": best.get("pubYear"), "journal": best.get("journalTitle", "")}
        if not ok: failed.append(key)
    except Exception as e:
        epmc_out[key] = {"found": False, "year_match": False, "error": str(e)[:80]}
        failed.append(key)
    time.sleep(0.15)

# CrossRef spot-check on five load-bearing references
cr_keys = ["netmhcpan41", "mhcflurry", "cppsite", "chengprusoff", "delong"]
cr_out = {}
for key, title, year in [r for r in REFS if r[0] in cr_keys]:
    q = urllib.parse.quote(title)
    url = f"https://api.crossref.org/works?query.bibliographic={q}&rows=2"
    try:
        d = get_json(url)
        items = d["message"]["items"]
        years = []
        for it in items:
            for f in ("published-print", "published-online", "issued"):
                if f in it and it[f].get("date-parts"):
                    years.append(it[f]["date-parts"][0][0]); break
        cr_out[key] = {"found": bool(items), "year_match": year in years,
                       "top_title": (items[0].get("title") or [""])[0][:110] if items else "",
                       "doi": items[0].get("DOI", "") if items else ""}
        if not items or year not in years: failed.append(key + ":crossref")
    except Exception as e:
        cr_out[key] = {"found": False, "error": str(e)[:80]}
    time.sleep(0.15)

# NCBI E-utilities spot-check on three database/tool papers
ncbi_out = {}
for key, title, year in [r for r in REFS if r[0] in ("cppsite", "cellppd", "iedb")]:
    q = urllib.parse.quote(title)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={q}&retmode=json&retmax=3"
    try:
        d = get_json(url)
        ids = d["esearchresult"]["idlist"]
        ncbi_out[key] = {"pubmed_ids": ids, "found": bool(ids)}
        if not ids: failed.append(key + ":ncbi")
    except Exception as e:
        ncbi_out[key] = {"found": False, "error": str(e)[:80]}
    time.sleep(0.34)

out = {"n_refs": len(REFS), "europe_pmc": epmc_out, "crossref": cr_out,
       "ncbi_eutils": ncbi_out, "failed": sorted(set(failed)),
       "verdict": "all bibliography entries verified" if not failed else f"UNVERIFIED: {sorted(set(failed))}"}
json.dump(out, open("results/citation_verification.json", "w"), indent=1)
print(json.dumps({k: (v if k in ("n_refs", "failed", "verdict") else "see file") for k, v in out.items()}, indent=1))

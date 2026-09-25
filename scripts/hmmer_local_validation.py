"""Reproduce 18 phmmer novelty searches against current full Swiss-Prot FASTA.

Usage: python scripts/hmmer_local_validation.py /path/to/phmmer /path/to/uniprot_sprot.fasta
The input FASTA is external and its provenance and checksum belong in the output.
"""
import hashlib, json, pathlib, subprocess, sys

binary, database = map(pathlib.Path, sys.argv[1:])
source = json.load(open('results/cpp_novel_candidates.json'))['named_candidates']
work = pathlib.Path('/tmp/09b-phmmer-local'); work.mkdir(exist_ok=True)
rows = {}
for c in source:
    name, sequence = c['name'], c['sequence']
    query = work / f'{name}.fa'; query.write_text(f'>{name}\n{sequence}\n')
    table, log = work / f'{name}.tbl', work / f'{name}.log'
    cmd = [str(binary), '--cpu', '2', '-E', '100', '--tblout', str(table),
           '--noali', '-o', str(log), str(query), str(database)]
    subprocess.run(cmd, check=True)
    lines = [x.split(maxsplit=18) for x in table.read_text().splitlines() if x and not x.startswith('#')]
    hits = [{'target': x[0], 'accession': x[1] if x[1] != '-' else x[0].split('|')[1], 'evalue': float(x[4]),
             'score': float(x[5]), 'description': x[18] if len(x)>18 else ''}
            for x in lines]
    rows[name] = {'sequence': sequence, 'n_hits_E_le_100': len(hits), 'top_hits': hits[:10]}
    print(name, 'reported<=100',len(hits),flush=True)
output = {'tool': 'HMMER 3.3.2 phmmer, Ubuntu 22.04 package 3.3.2+dfsg-1',
          'database_url': 'https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz',
          'database_download_date': '2026-09-25',
          'database_sha256_uncompressed': hashlib.sha256(database.read_bytes()).hexdigest(),
          'database_seq_count': sum(line.startswith(b'>') for line in database.open('rb')),
          'reporting_threshold_E': 100, 'results': rows}
pathlib.Path('results/cpp_hmmer_local_validation.json').write_text(json.dumps(output,indent=2)+'\n')

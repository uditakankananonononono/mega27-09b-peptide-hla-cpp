#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/raw
curl -sL --retry 3 -o data/raw/mhc_ligand_full_single_file.zip \
  "https://www.iedb.org/downloader.php?file_name=doc/mhc_ligand_full_single_file.zip"
curl -sL --retry 3 -o data/raw/cppsite2_natural.fa \
  "https://webs.iiitd.edu.in/raghava/cppsite/natural_pep.fa"
curl -sL --retry 3 -o data/raw/cppsite2_nonnatural.fa \
  "https://webs.iiitd.edu.in/raghava/cppsite/non-nat_pep.fa"
curl -sL --retry 3 -o data/raw/uniprot_reviewed_len8_35.fasta \
  "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28reviewed%3Atrue%29+AND+%28length%3A%5B8+TO+35%5D%29+AND+%28fragment%3Afalse%29"

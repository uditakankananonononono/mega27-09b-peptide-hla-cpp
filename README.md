# mega27-09b: Peptide-HLA binding prediction + cell-penetrating peptide design

MEGA-PROGRAM-27 item 9 (part 2).

## Projects
1. **pHLA-I binding prediction** on real IEDB MHC-ligand binding assays
   (quantitative IC50, class I, human alleles). CNN and GNN (residue-contact
   graph) models; benchmarked against published leader performance
   (NetMHCpan-4.1 / MHCflurry-2.0 reported metrics) on held-out alleles/peptides.
2. **Cell-penetrating peptide (CPP) classification + design** on real
   CPPsite 2.0 experimentally validated CPPs vs length-matched UniProt-derived
   negatives; classifier plus a conditional generator screened in silico
   (amphipathicity, net charge, hydrophobic moment, toxicity flags).

## Data sources (all real, open)
- IEDB MHC ligand export: https://www.iedb.org/downloader.php?file_name=doc/mhc_ligand_full_single_file.zip
- CPPsite 2.0 natural CPPs: https://webs.iiitd.edu.in/raghava/cppsite/natural_pep.fa
- UniProt reviewed fragments (negatives): https://rest.uniprot.org/uniprotkb/stream
- RCSB PDB peptide-HLA structures (graph anchors): https://files.rcsb.org/

## Rules (program-wide)
Real data only, no stubs/pseudocode. Hermetic pytest suite (live network only
outside CI). Math derivations + proofs in paper/. Honest benchmark verdicts.

## Command-line tool

The repository ships a working CLI over the trained, committed model
artifacts (`results/`). No training data download is needed.

```bash
pip install .                      # exposes the `peptidehlacpp` entry point
# or run in place: PYTHONPATH=src python3 -m peptidehlacpp.cli ...

# peptide-HLA class-I binding (PSSM log10 IC50 + CNN binder probability
# + the paper's z-scored ensemble for batches)
peptidehlacpp predict --allele 'HLA-A*02:01' --peptides GILGFVFTL AAAAAAAAA --json
peptidehlacpp predict --list-alleles

# CPP classification screen (trained CNN + biophysical descriptors)
peptidehlacpp screen --fasta candidates.fa --json

# novel CPP design: GRU generator + in-silico cascade
# (classifier >= 0.7, charge 2..12, hydrophobic moment, novelty vs train set)
peptidehlacpp design-cpp --n 10 --sample 1000 --seed 17 --json
```

Artifacts load from `results/` auto-detected upward from the package;
override with `--weights-dir`. Ensemble scores are z-scored within the
scored batch (as in the paper); single-peptide runs report PSSM/CNN only.
Tests: `python3 -m pytest tests/test_cli.py` (9 hermetic CLI tests).

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

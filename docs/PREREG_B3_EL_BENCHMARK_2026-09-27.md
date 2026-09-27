# PREREGISTRATION B3 - Eluted-ligand benchmark (locked 2026-09-27, BEFORE any
# B3 scoring)

Addendum item B3: "Ligand-eluted benchmark (NetMHCpan task) for a fair
comparison axis." G1's IEDB-tools-API egress block stands (403), so B3 uses
the offline NetMHCpan-4.1 training/evaluation data release as the locked
substitute axis.

## Data (source of record)
- File: NetMHCpan_train.tar.gz from
  https://services.healthtech.dtu.dk/suppl/immunology/NAR_NetMHCpan_NetMHCIIpan/NetMHCpan_train.tar.gz
  sha256 06f2c9f20bb959238bf5d601fca0489a0ed3f17648b952f30640205afca8f9b4
  (retrieved 2026-09-27; recorded in data/external/SOURCES.md; the tarball
  itself is NOT committed - 91 MB - but the hash + URL + retrieval date are).
- Cohort: UNION of the 5 EL partition files (c000_el..c004_el). We are not
  training on this data, so the CV partition structure is irrelevant; rows
  are deduplicated per (peptide, allele) with ligand (target 1) overriding
  decoy (target 0) on conflict. Locked.
- Mono-allelic rows only: rows whose 3rd column is a single HLA allele are
  IN; rows naming a multi-allelic cell line (allelelist-resolved) are
  EXCLUDED from the primary benchmark and counted + disclosed. Locked.
- Usable peptides: standard 20-letter alphabet, length 8-14.
- Allele name mapping: HLA-C14:02 -> HLA-C*14:02 (colon format to project
  format). Locked.

## Arms (both locked, both reported as-is)
- Primary (covered-allele arm): EL rows whose allele is in the 52-allele
  training set (results/per_allele_analysis.json). Scored by the FULL
  production ensemble (per-allele PSSM + allele-conditioned CNN, z-sum, the
  G1/G2 model family), reusing the B1-trained CNN (results/b1_partial/
  b1_cnn_best.pt, trained on all 52-allele IEDB BA rows - no EL row touches
  any fit).
- Secondary (uncovered-allele arm): EL rows whose allele is NOT in the
  training set, scored by the B1 arm-1 zero-shot mode (CNN + mean allele
  embedding). This quantifies EL-task transfer for unseen alleles, extending
  B1 from BA to EL data.

## Endpoints and decision rules (predeclared, B4-consistent)
- Per allele: AUROC (ligand vs decoy) and AUC0.1 (top-10% partial AUC, the
  NetMHCpan EL standard), allele-equal mean across alleles, 10,000-resample
  allele bootstrap, seed 23.
- The model is BA-trained; EL ranking is a transfer task. Locked bars
  (rationale disclosed here in advance): NetMHCpan-4.1 reports ~0.98 AUC on
  this task with EL training data; a BA-only model cannot be expected to
  match that, so -
  SUCCESS: mean AUROC >= 0.90 AND mean AUC0.1 >= 0.50.
  PARTIAL: 0.75 <= mean AUROC < 0.90.
  FAILURE: mean AUROC < 0.75.
- Single-class alleles (none expected - decoys are built in) would be
  excluded from means and disclosed, same locked handling as B1.
- Falsifier: seed-29 label shuffle within the largest covered allele must
  give AUROC in [0.40, 0.60]; otherwise the harness is void until fixed.
- Overlap audit: fraction of EL ligands that also appear in the IEDB BA
  training rows (peptide-level), reported per allele; flagged not removed.
  Locked.

## Reporting
results/b3_el_benchmark.json + per-allele table. Both arms as-is; the
outcome edits no prior locked result (G3, B1).

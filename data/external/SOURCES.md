# External data sources (pinned)
- MHC_pseudo.dat: MHC class I pseudo-sequences (34 contact positions per allele),
  as distributed with NetMHCpan (Nielsen et al., PLOS One 2007; NetMHCpan-4.1
  Reynisson et al., NAR 2020). Retrieved 2026-09-26 from the public mirror
  https://raw.githubusercontent.com/mnielLab/tcrcluster_backend/6c3449ca/data/Matrices/pseudoseqs_raw/MHC_pseudo.dat
  SHA256: see MHC_pseudo.sha256. Used for allele-held-out generalization (G3):
  allele features for alleles never seen in training.

## NetMHCpan-4.1 training release (B3)
- URL: https://services.healthtech.dtu.dk/suppl/immunology/NAR_NetMHCpan_NetMHCIIpan/NetMHCpan_train.tar.gz
- sha256: 06f2c9f20bb959238bf5d601fca0489a0ed3f17648b952f30640205afca8f9b4
- retrieved 2026-09-27; tarball + extracted files NOT committed (size); regenerable from URL.

# Paper-build evidence log, 2026-09-27

Branch: `paper-build`, paper directory only. `paper/manuscript.md` adapts the existing root `paper_full_draft.md` without changing that legacy file.

- Corrected an invalid identifiability theorem: full marginal amino-acid coverage and a large sample size do not imply full rank of a peptide one-hot design. The corrected statement requires actual full rank after choosing a gauge; an empirical rank audit remains pending. Removed a claim that the observed design covariance must be block diagonal.
- Corrected a contradictory CNN AUC citation (0.9057 to the draft's result-table value 0.9044). Added the negative locked G3 allele-held-out result from `results/g3_allele_holdout.json`: AUROC 0.6687 versus MHCflurry 0.9285 on 6,000 examples, G3 gate failed.
- Compiled `paper/manuscript.pdf` from `manuscript.md` via `pandoc -t latex` and pdfLaTeX: six pages by `pdfinfo`. It is **not** a 50-page submission-ready manuscript and is thin in methods, figures, and references. The body uses a Times-compatible Nimbus Roman font (`mathptmx`), not genuine TNR. Visual check: readable but equations inherited from the draft are plain-text approximations, so typesetting remains open. Do not promote it as a final paper.

## Follow-up expansion, 2026-09-27

- Added `companion_analysis.md` (rendered into `companion_body.tex`): same-row MHCflurry table with paired bootstrap, G2 per-allele comparison, all twelve G3 held-out allele AUCs, three assay-method slices, confound/calibration caveats, and 18 computational CPP nomination rows, each from committed result JSON. The within-panel numerical win is explicitly subordinate to G3's transfer failure.
- Rebuilt `manuscript.pdf` with pdfLaTeX twice: 10 pages, up from six. `mathptmx` resolves to Nimbus Roman in this environment, not genuine TNR. Visual inspection of pages 8-10 found tables readable after shortening labels and separating long sequences from metrics. Formula prose from the original draft is still not publication-grade; 50-page, full reference and font gates remain open.

## Morning evidence expansion, 2026-09-27

- Added `detailed_diagnostics.md` (compiled into `detailed_body.tex`) from committed per-allele head-to-head, per-length head-to-head, CPP learning curve and k-mer feature-weight artifacts. The 30 familiar-allele rows and length slices are deliberately separate from the failed 12-allele G3 transfer test; small strata cannot support a claim of mechanism.
- Rebuilt with pdfLaTeX twice: 13 pages, up from ten, Nimbus Roman substitute rather than licensed Times New Roman. Visual inspection of pages 11-13 showed readable metrics and no material table clipping. Still not a 50-page paper; do not add blank or repetitive pages to meet a numerical floor.

# Paper-build evidence log, 2026-09-27

Branch: `paper-build`, paper directory only. `paper/manuscript.md` adapts the existing root `paper_full_draft.md` without changing that legacy file.

- Corrected an invalid identifiability theorem: full marginal amino-acid coverage and a large sample size do not imply full rank of a peptide one-hot design. The corrected statement requires actual full rank after choosing a gauge; an empirical rank audit remains pending. Removed a claim that the observed design covariance must be block diagonal.
- Corrected a contradictory CNN AUC citation (0.9057 to the draft's result-table value 0.9044). Added the negative locked G3 allele-held-out result from `results/g3_allele_holdout.json`: AUROC 0.6687 versus MHCflurry 0.9285 on 6,000 examples, G3 gate failed.
- Compiled `paper/manuscript.pdf` from `manuscript.md` via `pandoc -t latex` and pdfLaTeX: six pages by `pdfinfo`. It is **not** a 50-page submission-ready manuscript and is thin in methods, figures, and references. The body uses a Times-compatible Nimbus Roman font (`mathptmx`), not genuine TNR. Visual check: readable but equations inherited from the draft are plain-text approximations, so typesetting remains open. Do not promote it as a final paper.

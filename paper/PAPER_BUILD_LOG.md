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

## Structure and error audit, 2026-09-27

- Added `structure_audit.tex` from committed HLA pocket-entropy (52 alleles, 23 sites), MDAnalysis groove-geometry (69 analyzed complexes, 3 excluded), three selected anchor-burial structures, motif interaction summaries, PSSM learning curve and paired-error counts. The long structure table is an audit of distinct deposited complexes, not an independent binding replication. The pocket-residue-66 group comparison is negative (Kruskal p=0.20685), while the small within-panel pair-error edge cannot overturn G3.
- Rebuilt twice after removing an accidental duplicate include: 18 pages, up from 13. Visually inspected pages 14-18; tables are readable and not clipped. Page 18 is sparse because the previous longtable and the final short diagnostic paragraph cross a page break, not because blank pages were inserted. Nimbus Roman, not genuine Times New Roman; the 50-page floor remains unmet.

## Judge requirement amended, 2026-09-27 10:00 IST

The owner said "NOT 10 ROUNDS OOF CHATGPT CHECK JUST ONE WHICH I PROVIDE OK?" (authenticated WhatsApp message `wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMDJCMTZGRTVEMkQwMTFBQzc4MQA=`, 10:00:07 IST). For this project, the paper branch therefore marks the counted judge gate **0 of 1, PENDING her personally provided verdict**. A round initiated by agents, even through her ChatGPT account, remains historical or supplementary and does not meet the gate. Historical ten-round language in the science ledgers is not erased by this note. A project-specific user-pasted verdict must be traced and evaluated before completion is recorded. Supplementary Gemini/LLM consults do not count. No scientific result, page count or font gate changes here.

## Authorship-attribution cleanup, 2026-09-27 11:14 IST

The owner requested removal of the assistant's attribution from the papers (WhatsApp `wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMEY5MzY4M0Q4OUYwNjg4ODZDNwA=`). Removed agent/program-style byline and credit text from the editable paper source and PDF display, without substituting an author. Udita's own byline in 09b was preserved, with only the Instinct pipeline parenthetical removed. Manuscript PDF author metadata is empty. Literature references to other studies' authors and technical uses of "author numbering" are not authorship credits for this paper.

## Owner mega-verdict, 2026-09-27 noon IST

Original complete WhatsApp message archived in `MEGA_VERDICT_FULL_BODY_2026-09-27.txt` (SHA-256 `d700f12c2a01d6f21b7392305aaa6b25d7a9efb29062ce5ba95c14e22f11bc33`); lane and cross-cutting excerpts in `MEGA_VERDICT_2026-09-27.md`, with a separate agent-authored locked queue. Source: wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMEFGRDY4MzY4OTkxNzFEQURGRAA=.

Added one-question framing and a G3-first interpretive response. Committed G3 result remains negative: 0.6687 against MHCflurry 0.9285 pooled, and 0.5798 versus 0.9056 mean per allele. Mechanisms are hypotheses, not tested explanations; no new benchmark numbers invented. PDF double-built to 20 pages, with Nimbus Roman substitute, not licensed Times New Roman; visual inspection of final two pages completed. The 50-page font gate remains open. No slides rendered in this paper-side revision.

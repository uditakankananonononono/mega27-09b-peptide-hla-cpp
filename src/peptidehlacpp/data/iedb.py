"""Streaming filter for the IEDB MHC-ligand full export.

The raw export is a 9.2GB single-file CSV inside a zip. We never materialize
it: rows are streamed, two header rows are parsed, and only class I human
HLA rows with a quantitative nM affinity and a standard 8-15mer peptide are
kept. Output is a compact TSV used by the modeling pipeline.
"""
from __future__ import annotations

import csv
import io
import re
import sys
import zipfile
from pathlib import Path

AA = set("ACDEFGHIKLMNPQRSTVWY")

# Column indices in mhc_ligand_full.csv (verified against the 2026-09-22 export)
COL_EPITOPE_SEQ = 11      # Epitope / Name
COL_HOST_NAME = 43        # Host / Name
COL_ASSAY_METHOD = 90     # Assay / Method
COL_ASSAY_UNITS = 92      # Assay / Units
COL_QUALITATIVE = 94      # Assay / Qualitative Measurement
COL_INEQUALITY = 95       # Assay / Measurement Inequality
COL_QUANTITATIVE = 96     # Assay / Quantitative measurement
COL_MHC_NAME = 107        # MHC Restriction / Name
COL_MHC_CLASS = 111       # MHC Restriction / Class
COL_REFERENCE_IRI = 1     # Reference / IEDB IRI (row groups share names; idx 1 is Reference IRI)

HLA_RE = re.compile(r"HLA-[ABC]\*\d{2}:\d{2}")


def normalize_allele(raw: str) -> str | None:
    """Extract a canonical HLA-?*NN:NN allele name, or None if not class-I HLA."""
    m = HLA_RE.search(raw or "")
    return m.group(0) if m else None


def row_is_usable(fields: list[str]) -> tuple[bool, str, str, str, float, str, str]:
    """Return (keep, sequence, allele, units, value, qualitative, reference)."""
    seq = fields[COL_EPITOPE_SEQ].strip().upper()
    if not (8 <= len(seq) <= 15) or not set(seq) <= AA:
        return (False, "", "", "", 0.0, "", "")
    allele = normalize_allele(fields[COL_MHC_NAME])
    if allele is None or fields[COL_MHC_CLASS].strip() != "I":
        return (False, "", "", "", 0.0, "", "")
    units = fields[COL_ASSAY_UNITS].strip()
    if units != "nM":
        return (False, "", "", "", 0.0, "", "")
    try:
        value = float(fields[COL_QUANTITATIVE])
    except (ValueError, IndexError):
        return (False, "", "", "", 0.0, "", "")
    if not (0.0 < value <= 1e7):
        return (False, "", "", "", 0.0, "", "")
    qual = fields[COL_QUALITATIVE].strip()
    ref = fields[COL_REFERENCE_IRI].strip()
    return (True, seq, allele, units, value, qual, ref)


def filter_export(zip_path: str | Path, out_tsv: str | Path,
                  log_every: int = 2_000_000) -> dict:
    """Stream-filter the IEDB export into a compact TSV. Returns stats."""
    stats = {"rows_scanned": 0, "rows_kept": 0}
    out_fields = ["sequence", "allele", "units", "value_nm",
                  "inequality", "qualitative", "method", "host", "reference_iri"]
    with zipfile.ZipFile(zip_path) as zf, \
         zf.open("mhc_ligand_full.csv") as raw, \
         io.TextIOWrapper(raw, encoding="utf-8", errors="replace") as text, \
         open(out_tsv, "w") as out:
        reader = csv.reader(text)
        next(reader)  # group header row
        next(reader)  # column header row
        out.write("\t".join(out_fields) + "\n")
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        for fields in reader:
            stats["rows_scanned"] += 1
            if len(fields) <= COL_MHC_CLASS:
                continue
            keep, seq, allele, units, value, qual, ref = row_is_usable(fields)
            if keep:
                stats["rows_kept"] += 1
                writer.writerow([seq, allele, units, value,
                                 fields[COL_INEQUALITY].strip(), qual,
                                 fields[COL_ASSAY_METHOD].strip(),
                                 fields[COL_HOST_NAME].strip(), ref])
            if log_every and stats["rows_scanned"] % log_every == 0:
                print(f"scanned={stats['rows_scanned']} kept={stats['rows_kept']}",
                      file=sys.stderr, flush=True)
    return stats


if __name__ == "__main__":
    import json
    zp = sys.argv[1] if len(sys.argv) > 1 else "data/raw/mhc_ligand_full_single_file.zip"
    op = sys.argv[2] if len(sys.argv) > 2 else "data/processed/iedb_class1_human_nM.tsv"
    Path(op).parent.mkdir(parents=True, exist_ok=True)
    print(json.dumps(filter_export(zp, op)))

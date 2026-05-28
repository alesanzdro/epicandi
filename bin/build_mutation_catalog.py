#!/usr/bin/env python3
"""
build_mutation_catalog.py — derive a clean, panel-aware mutation catalog
from the FungAMR CSV (filtered to a given species / host).

The raw FungAMR CSV (~27 500 rows, 30 cols) carries:
  - multi-gene rows: 'gene or protein name' may be 'Cdr1,Erg11,Fcy1,Fks1'
    meaning the mutation was reported in a strain carrying changes in all
    of those genes — we expand to one row per gene.
  - heterogeneous mutation strings: 'V125A' / 'F635del' / 'E86fs' / 'Y177X'
    / 'Knockout' / 'Chr5x2' / 'Translocation'. We classify each and parse
    SNP / Indel / LoF when the syntax allows.
  - confidence_score: positive = evidence FOR resistance (1=strongest,
    8=weakest), negative = evidence AGAINST resistance. We summarise per
    (gene, mutation) as the most extreme score seen across reports.

The output is a flat catalog plugged into:
  - aggregate_cnv.py (R/I/S modulation),
  - cnv_amr_atlas.py (annotate SNPs / CNV calls with FungAMR evidence),
  - extract_gene_proteins.py (classify per-gene mismatches).

Output: TSV with columns
  gene, mutation_id, mutation_type, ref_aa, position, mut_aa,
  drugs, confidence_score, evidence_strength, n_reports, pmids,
  tier, in_panel

Filter knobs:
  --species  "auris"   substring match against FungAMR 'species' col (case-insensitive)
  --host     ""        substring match against FungAMR 'host' col (empty = no filter)
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


SNP_RE = re.compile(r"^([A-Z*])\s*([0-9]+)\s*([A-Z*])$")
LOF_RE = re.compile(r"^([A-Z])\s*([0-9]+)\s*(X|\*)$")             # Y177X / Y177*
FS_RE = re.compile(r"^([A-Z])\s*([0-9]+)\s*(fs|FS|frameshift)$")  # E86fs
DEL_RE = re.compile(r"^([A-Z]?)\s*([0-9]+)\s*(del|Δ|DEL)$", re.I) # F635del / Δ33
RANGE_DEL_RE = re.compile(r"^([A-Z])([0-9]+)_([A-Z])([0-9]+)(del|Δ|DEL)$", re.I)  # M1_K33del


def split_compound_mutation(mut: str) -> list[str]:
    """A FungAMR row may report multiple mutations seen jointly in one isolate
    (e.g. 'V704L|K143R|S70R'). Split on common separators."""
    if not mut:
        return []
    # Don't split 'M1_K33del' style ranges (underscore is part of the token).
    parts = re.split(r"[|,;]\s*", mut)
    return [p.strip() for p in parts if p.strip()]


def classify_mutation(mut_str: str, mtype_raw: str) -> tuple[str, str, int | None, str]:
    """Return (mutation_type_clean, ref_aa, position, mut_aa)."""
    mut = (mut_str or "").strip()
    raw = (mtype_raw or "").strip().lower()

    # Highest-fidelity: parse the syntax of the mutation string.
    if m := SNP_RE.match(mut):
        return "SNP", m.group(1), int(m.group(2)), m.group(3)
    if m := LOF_RE.match(mut):
        return "LoF", m.group(1), int(m.group(2)), "X"
    if m := FS_RE.match(mut):
        return "Indel", m.group(1), int(m.group(2)), "fs"
    if m := RANGE_DEL_RE.match(mut):
        # M1_K33del → record the first position as the anchor
        return "Indel", m.group(1), int(m.group(2)), f"del_{m.group(3)}{m.group(4)}"
    if m := DEL_RE.match(mut):
        ref = m.group(1) or ""
        return "Indel", ref, int(m.group(2)), "del"

    # Fall back to the FungAMR-supplied type label.
    if "cnv" in raw or "knockout" in mut.lower() or "amplification" in mut.lower() or "duplication" in mut.lower():
        return "CNV", "", None, mut
    if "loss" in raw or "function" in raw:
        return "LoF", "", None, mut
    if "indel" in raw:
        return "Indel", "", None, mut
    if "snp" in raw:
        return "SNP", "", None, mut
    return "Other", "", None, mut


def evidence_label(score: float | None) -> str:
    if score is None:
        return "NA"
    if score > 0:
        return "R-strong" if score <= 4 else "R-weak"
    if score < 0:
        return "S-strong" if score >= -4 else "S-weak"
    return "neutral"


def normalize_drug(d: str) -> str:
    s = (d or "").strip().lower()
    if not s:
        return ""
    # 5-Flucytosine / 5FC / Flucytosine variants
    if "flucytosine" in s or s == "5fc":
        return "5fc"
    if s == "amphotericin b" or s == "amphotericin_b" or s == "ampb":
        return "amphotericin_b"
    if "echinocandin" in s:
        return "echinocandins"
    # else: lowercase, spaces → underscores
    return s.replace(" ", "_")


def load_panel_tiers(panel_path: str | None) -> dict[str, str]:
    """gene → tier; empty dict if no panel."""
    if not panel_path or not Path(panel_path).exists():
        return {}
    out = {}
    with open(panel_path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            g = (row.get("gene") or "").strip().upper()
            if g:
                out[g] = (row.get("tier") or "").strip()
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fungamr", required=True,
                    help="FungAMR CSV (e.g. resources/databases/fungamr/FungAMR_070425_epicandi.csv)")
    ap.add_argument("--species", default="auris",
                    help="Case-insensitive substring of FungAMR 'species' column (default: auris)")
    ap.add_argument("--host", default="",
                    help="Case-insensitive substring of FungAMR 'host' column (default: no filter)")
    ap.add_argument("--panel", default=None,
                    help="assets/cnv_loci/panel.tsv (for tier annotation)")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    panel_tiers = load_panel_tiers(args.panel)
    species_q = args.species.lower()
    host_q = args.host.lower()

    # group key = (gene_upper, mutation_string_normalized)
    # aggregate: drugs (set), pmids (set), scores (list), n_reports,
    # companion_mutations (set of mutations seen jointly in the same isolate)
    grouped: dict[tuple[str, str], dict] = defaultdict(lambda: {
        "drugs": set(), "pmids": set(), "scores": [], "n_reports": 0,
        "mtype_raw": "", "mut_orig": "", "companions": set(),
    })

    n_total = n_kept = n_expanded = 0
    with open(args.fungamr, encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            n_total += 1
            sp = (row.get("species") or "").lower()
            host = (row.get("host") or "").lower()
            if species_q and species_q not in sp:
                continue
            if host_q and host_q not in host:
                continue
            n_kept += 1

            gene_field = (row.get("gene or protein name") or "").strip()
            mut = (row.get("mutation") or "").strip()
            drug = normalize_drug(row.get("drug", ""))
            try:
                score = float(row.get("confidence score") or "")
            except ValueError:
                score = None
            pmid = (row.get("pubmedid") or "").strip()
            mtype_raw = (row.get("mutation_type") or "").strip()

            # Compound mutation strings (joint isolate) → split AND keep the
            # rest as companions so we can flag 'this SNP is reported as causal
            # only when accompanied by X|Y|Z'.
            sub_muts = split_compound_mutation(mut) or [""]

            for g_raw in gene_field.split(","):
                g = g_raw.strip().upper()
                if not g:
                    continue
                for sub_mut in sub_muts:
                    n_expanded += 1
                    key = (g, sub_mut)
                    rec = grouped[key]
                    rec["n_reports"] += 1
                    if drug:
                        rec["drugs"].add(drug)
                    if pmid:
                        rec["pmids"].add(pmid)
                    if score is not None:
                        rec["scores"].append(score)
                    if not rec["mtype_raw"]:
                        rec["mtype_raw"] = mtype_raw
                    if not rec["mut_orig"]:
                        rec["mut_orig"] = mut
                    # companions = the other mutations in the same isolate
                    companions = [m for m in sub_muts if m and m != sub_mut]
                    for c in companions:
                        rec["companions"].add(c)

    print(f"[catalog] FungAMR rows total: {n_total}",                file=sys.stderr)
    print(f"[catalog] after filter species='{args.species}' host='{args.host}': {n_kept}",
          file=sys.stderr)
    print(f"[catalog] after multi-gene expansion: {n_expanded}",     file=sys.stderr)
    print(f"[catalog] unique (gene, mutation) groups: {len(grouped)}",file=sys.stderr)

    cols = ["gene", "mutation_id", "mutation_type", "ref_aa", "position",
            "mut_aa", "drugs", "confidence_score", "evidence_strength",
            "n_reports", "pmids", "companion_mutations", "tier", "in_panel"]
    out_rows = []
    for (g, mut), rec in grouped.items():
        mtype, ref_aa, pos, mut_aa = classify_mutation(mut, rec["mtype_raw"])
        # most-extreme positive score (preferred); else most-extreme negative
        scores_pos = [s for s in rec["scores"] if s > 0]
        scores_neg = [s for s in rec["scores"] if s < 0]
        if scores_pos:
            score_repr = min(scores_pos)   # lower = stronger evidence for R
        elif scores_neg:
            score_repr = max(scores_neg)   # closer to 0 = stronger evidence for R (weaker S)
        else:
            score_repr = None
        out_rows.append({
            "gene": g,
            "mutation_id": mut or rec["mut_orig"],
            "mutation_type": mtype,
            "ref_aa": ref_aa,
            "position": pos if pos is not None else "",
            "mut_aa": mut_aa,
            "drugs": ";".join(sorted(rec["drugs"])),
            "confidence_score": f"{score_repr}" if score_repr is not None else "",
            "evidence_strength": evidence_label(score_repr),
            "n_reports": rec["n_reports"],
            "pmids": ";".join(sorted(rec["pmids"])),
            "companion_mutations": ";".join(sorted(rec["companions"])),
            "tier": panel_tiers.get(g, ""),
            "in_panel": "true" if g in panel_tiers else "false",
        })

    # Sort by gene, then by position (numeric where available), then mutation_id
    def sort_key(r):
        pos = r["position"] if r["position"] != "" else 1e9
        return (r["gene"], pos, r["mutation_id"])
    out_rows.sort(key=sort_key)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in out_rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"[catalog] wrote {args.output} ({len(out_rows)} mutations)", file=sys.stderr)

    # Quick stats
    by_type = defaultdict(int)
    in_panel_total = 0
    for r in out_rows:
        by_type[r["mutation_type"]] += 1
        if r["in_panel"] == "true":
            in_panel_total += 1
    print(f"[catalog] by type: " + ", ".join(f"{k}={v}" for k, v in sorted(by_type.items())),
          file=sys.stderr)
    print(f"[catalog] in panel: {in_panel_total}/{len(out_rows)}", file=sys.stderr)


if __name__ == "__main__":
    main()

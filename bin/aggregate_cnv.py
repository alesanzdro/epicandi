#!/usr/bin/env python3
"""
aggregate_cnv.py — cohort-wide CNV aggregator with mosdepth-CNVkit consensus
and optional AMR modulation by CNV events.

Inputs (globs accept multiple paths via shell expansion or repeated flags):
  --events       *.cnv_events.tsv      (one per sample, from cnv_detect.py)
  --cnvkit_panel *.cnvkit_panel.tsv    (one per sample, from cnvkit_extract_panel.py)
  --contig_depth *.contig_depth.tsv    (one per sample, from cnv_detect.py)
  --call_matrix  call_matrix.tsv       (cohort-wide AMR from ChroQueTas, optional)
  --cnv_amr_rules cnv_amr_rules.tsv    (assets/cnv_loci/, only used if --call_matrix)
  --panel        panel.tsv             (gene order, tiers)
  --outdir       results/cnv/aggregated/

Outputs in --outdir:
  cnv_call_matrix_mosdepth.tsv      samples x genes (event labels)
  cnv_call_matrix_cnvkit.tsv        samples x genes (event labels)
  cnv_call_matrix_consensus.tsv     samples x genes (consensus or LOW_CONFIDENCE)
  cnv_log2_matrix.tsv               samples x genes (mosdepth log2 ratio)
  cnv_events_summary.tsv            long format, one row per non-normal event
  cnv_contig_depth_combined.tsv     long format, all contigs all samples
  call_matrix_with_cnv.tsv          AMR call matrix modulated by CNV
                                    (only if --call_matrix is provided)

Consensus logic per (sample, gene):
  Both callers agree at the bucket level (dup/del/normal/no_cov)
      -> consensus = the mosdepth event (preserves granularity:
         amplification > dup_high > dup_low)
  Buckets disagree
      -> consensus = LOW_CONFIDENCE_<mosdepth_event>_<cnvkit_event>
  Only one caller has coverage
      -> consensus = <event>__single_caller
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


# Mapping mosdepth events to CNVkit's coarser bucket
MOSDEPTH_TO_BUCKET = {
    "amplification": "dup",
    "dup_high": "dup",
    "dup_low": "dup",
    "normal": "normal",
    "deletion": "del",
    "deletion_total": "del",
    "no_coverage": "no_cov",
    "aneuploidy_2x": "dup",
    "aneuploidy_3x": "dup",
}


def read_tsv_rows(path: str) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_panel_gene_order(path: str) -> list[str]:
    out: list[str] = []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            g = row.get("gene", "").strip()
            if g and not g.startswith("#"):
                out.append(g)
    return out


def load_rules(path: str) -> list[dict]:
    if not path or not Path(path).exists():
        return []
    out = []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if not row.get("gene") or row["gene"].startswith("#"):
                continue
            out.append(row)
    return out


def write_matrix(path: str, samples: list[str], genes: list[str],
                 data: dict[tuple[str, str], str], default: str = "") -> None:
    with open(path, "w") as fh:
        fh.write("sample_id\t" + "\t".join(genes) + "\n")
        for s in samples:
            row = [s] + [data.get((s, g), default) for g in genes]
            fh.write("\t".join(row) + "\n")


def consensus_call(mos_event: str, cnv_event: str) -> str:
    if mos_event == "no_coverage" and cnv_event == "no_cov":
        return "no_coverage"
    if mos_event == "no_coverage":
        return f"{cnv_event}__single_caller"
    if cnv_event == "no_cov":
        return f"{mos_event}__single_caller"
    mb = MOSDEPTH_TO_BUCKET.get(mos_event, "normal")
    cb = cnv_event  # already a bucket
    if mb == cb:
        return mos_event  # keep the more granular mosdepth label
    return f"LOW_CONFIDENCE_{mos_event}_{cnv_event}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", nargs="+", required=True,
                    help="One or more *.cnv_events.tsv files")
    ap.add_argument("--cnvkit_panel", nargs="*", default=[],
                    help="One or more *.cnvkit_panel.tsv files (optional)")
    ap.add_argument("--contig_depth", nargs="*", default=[],
                    help="One or more *.contig_depth.tsv files")
    ap.add_argument("--call_matrix", default=None,
                    help="AMR call_matrix.tsv from ChroQueTas (optional)")
    ap.add_argument("--cnv_amr_rules", default=None,
                    help="cnv_amr_rules.tsv (only used if --call_matrix is given)")
    ap.add_argument("--panel", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    genes = load_panel_gene_order(args.panel)

    # ── load mosdepth events ──
    mos_event: dict[tuple[str, str], str] = {}
    mos_log2: dict[tuple[str, str], str] = {}
    samples: set[str] = set()
    summary_rows: list[dict] = []
    for path in args.events:
        for row in read_tsv_rows(path):
            sid, gene = row["sample_id"], row["gene"]
            samples.add(sid)
            mos_event[(sid, gene)] = row["event"]
            mos_log2[(sid, gene)] = row.get("log2_ratio", "")
            if row["event"] not in ("normal", "no_coverage"):
                summary_rows.append({
                    "sample_id": sid, "gene": gene, "tier": row.get("tier", ""),
                    "event": row["event"], "norm": row.get("norm", ""),
                    "log2_ratio": row.get("log2_ratio", ""),
                    "confidence": row.get("confidence", ""),
                    "expected_event": row.get("expected_event", ""),
                    "agreement_with_expected": row.get("agreement_with_expected", ""),
                    "source": "mosdepth",
                })

    # ── load CNVkit events ──
    cnv_event: dict[tuple[str, str], str] = {}
    for path in args.cnvkit_panel:
        for row in read_tsv_rows(path):
            sid, gene = row["sample_id"], row["gene"]
            samples.add(sid)
            cnv_event[(sid, gene)] = row["cnvkit_event"]
            if row["cnvkit_event"] not in ("normal", "no_cov"):
                summary_rows.append({
                    "sample_id": sid, "gene": gene, "tier": "",
                    "event": row["cnvkit_event"],
                    "norm": "", "log2_ratio": row.get("cnvkit_log2", ""),
                    "confidence": "", "expected_event": "",
                    "agreement_with_expected": "",
                    "source": "cnvkit",
                })

    samples_sorted = sorted(samples)
    panel_genes = [g for g in genes if any((s, g) in mos_event or (s, g) in cnv_event
                                           for s in samples_sorted)]

    # ── write matrices ──
    write_matrix(str(outdir / "cnv_call_matrix_mosdepth.tsv"),
                 samples_sorted, panel_genes, mos_event, default="no_coverage")
    write_matrix(str(outdir / "cnv_call_matrix_cnvkit.tsv"),
                 samples_sorted, panel_genes, cnv_event, default="no_cov")
    write_matrix(str(outdir / "cnv_log2_matrix.tsv"),
                 samples_sorted, panel_genes, mos_log2, default="")

    consensus: dict[tuple[str, str], str] = {}
    n_low_conf = 0
    for s in samples_sorted:
        for g in panel_genes:
            mev = mos_event.get((s, g), "no_coverage")
            cev = cnv_event.get((s, g), "no_cov")
            c = consensus_call(mev, cev)
            consensus[(s, g)] = c
            if c.startswith("LOW_CONFIDENCE_"):
                n_low_conf += 1
    write_matrix(str(outdir / "cnv_call_matrix_consensus.tsv"),
                 samples_sorted, panel_genes, consensus, default="no_coverage")

    # ── events summary (long) ──
    summary_path = outdir / "cnv_events_summary.tsv"
    with open(summary_path, "w") as fh:
        cols = ["sample_id", "gene", "tier", "event", "norm", "log2_ratio",
                "confidence", "expected_event", "agreement_with_expected", "source"]
        fh.write("\t".join(cols) + "\n")
        for r in summary_rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    # ── contig depth combined ──
    contig_path = outdir / "cnv_contig_depth_combined.tsv"
    cols_contig = ["sample_id", "contig", "length_bp", "mean_depth",
                   "ratio_vs_genome", "log2_ratio", "event"]
    with open(contig_path, "w") as fh:
        fh.write("\t".join(cols_contig) + "\n")
        for path in args.contig_depth:
            with open(path) as ifh:
                next(ifh, None)  # skip header
                for line in ifh:
                    fh.write(line if line.endswith("\n") else line + "\n")

    print(f"[aggregate_cnv] {len(samples_sorted)} samples x {len(panel_genes)} genes",
          file=sys.stderr)
    print(f"[aggregate_cnv] LOW_CONFIDENCE cells: {n_low_conf}", file=sys.stderr)
    print(f"[aggregate_cnv] non-normal events (long format): {len(summary_rows)}",
          file=sys.stderr)
    print(f"[aggregate_cnv] outputs in {outdir}/", file=sys.stderr)

    # ── AMR ↔ CNV modulation ──
    if args.call_matrix and Path(args.call_matrix).exists():
        rules = load_rules(args.cnv_amr_rules or "")
        if not rules:
            print(f"[aggregate_cnv] no rules loaded — skipping AMR modulation",
                  file=sys.stderr)
            return

        with open(args.call_matrix) as fh:
            amr_header = fh.readline().rstrip("\n").split("\t")
            amr_rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]

        # The AMR matrix layout: first column = sample_id, the rest = drug or
        # drug__call columns. We treat any column whose name appears in the
        # rules `drug` field as modulatable.
        amr_drug_cols = {c for c in amr_header[1:] if any(c == r["drug"] for r in rules)}
        modulation_log: list[dict] = []

        amr_out_rows = []
        amr_out_header = amr_header + ["modulation_source", "modulated_drugs",
                                       "amr_call_snp_only"]
        for row in amr_rows:
            sid = row[0]
            row_map = dict(zip(amr_header, row))
            sources = []
            modulated = []
            snp_only_snapshot = {d: row_map.get(d, "") for d in amr_drug_cols}

            for rule in rules:
                gene = rule["gene"]
                event_req = rule["event"]
                drug = rule["drug"]
                mod = rule["modulation"]
                mos_call = mos_event.get((sid, gene), "no_coverage")
                consensus_call_value = consensus.get((sid, gene), "no_coverage")
                # Trigger if either mosdepth call or consensus matches; require
                # that consensus is NOT LOW_CONFIDENCE_ unless mosdepth call
                # alone is strong evidence.
                triggered = (mos_call == event_req or consensus_call_value == event_req)
                if not triggered:
                    continue
                if drug not in row_map:
                    continue
                current = row_map[drug]
                new_value = current
                if mod == "override_to_S":
                    new_value = "S"
                elif mod == "override_to_R":
                    new_value = "R"
                elif mod == "reinforce_R":
                    if current in ("S", "I", ""):
                        new_value = "R"
                    # else: keep current=R, still register as supporting evidence
                elif mod == "flag_for_review":
                    new_value = f"{current}*" if current else "?*"

                # Always record the rule's application (even if value unchanged)
                # so clinicians can trace CNV-derived support for the call.
                row_map[drug] = new_value
                sources.append(f"{gene}:{event_req}:{mod}")
                if new_value != current:
                    modulated.append(drug)
                modulation_log.append({
                    "sample_id": sid, "drug": drug,
                    "rule": f"{gene}:{event_req}:{mod}",
                    "before": current, "after": new_value,
                    "source": rule.get("source", ""),
                })

            row_map["modulation_source"] = ";".join(sources)
            row_map["modulated_drugs"] = ";".join(modulated)
            row_map["amr_call_snp_only"] = ";".join(
                f"{d}={v}" for d, v in snp_only_snapshot.items() if v)
            amr_out_rows.append([row_map.get(c, "") for c in amr_out_header])

        amr_out_path = outdir / "call_matrix_with_cnv.tsv"
        with open(amr_out_path, "w") as fh:
            fh.write("\t".join(amr_out_header) + "\n")
            for r in amr_out_rows:
                fh.write("\t".join(r) + "\n")

        log_path = outdir / "cnv_amr_modulation_log.tsv"
        with open(log_path, "w") as fh:
            cols_log = ["sample_id", "drug", "rule", "before", "after", "source"]
            fh.write("\t".join(cols_log) + "\n")
            for r in modulation_log:
                fh.write("\t".join(str(r.get(c, "")) for c in cols_log) + "\n")

        print(f"[aggregate_cnv] AMR↔CNV modulations: {len(modulation_log)}",
              file=sys.stderr)
        print(f"[aggregate_cnv] wrote {amr_out_path} and {log_path}", file=sys.stderr)


if __name__ == "__main__":
    main()

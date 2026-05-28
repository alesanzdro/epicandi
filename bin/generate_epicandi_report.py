#!/usr/bin/env python3
"""
generate_epicandi_report.py
===========================

Builds the EPIMOL/FISABIO HTML run report for an epicandi cohort.

Phases implemented:
  0 + 1.  Banner + summary cards (samples, GB, PASS, AMR, CNV)
  2.      Per-species batch composition + identification chart
  3.      Fastq QC section (raw vs clean, 4-panel ANI, FAIL table) + secondary HTML
  4.      BUSCO + QUAST assembly metrics
  5.      AMR + CNV heatmaps (embedding existing PNGs)
  6.      SNP distance heatmap + UPGMA network (embed SVG)
  7.      Excel companion (6 sheets)

Spec: ~/.claude/projects/-home-asanzc-epicandi/memory/project_html_report_spec.md

Usage
-----
    python3 bin/generate_epicandi_report.py \\
        --results-dir results_hospital60 \\
        --run-name 260526_HGRAL_HLAFE \\
        --samplesheet test/samplesheet_hospital60.csv \\
        --branding branding/ \\
        --output results_hospital60/00_report/260526_HGRAL_HLAFE_epicandi_report.html
"""
from __future__ import annotations

import argparse
import base64
import csv
import io
import json
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np


# ─── PALETTE ─────────────────────────────────────────────────────────────────
BLUE   = "#2C5F8A"
LBLUE  = "#5B9ABF"
ORANGE = "#E05C2A"
GREEN  = "#28A745"
RED    = "#E74C3C"
GREY   = "#6C757D"
DARK   = "#1F2933"
BG     = "#F4F6F9"

# Cycle used to colour distinct batches. Extended beyond the corporate palette
# so cohorts with many batches still get visually distinct colours.
BATCH_PALETTE = [
    BLUE, ORANGE, GREEN, "#9B59B6", "#F39C12",
    "#1ABC9C", RED, "#3498DB", "#95A5A6", "#D35400",
    "#16A085", "#2980B9", "#C0392B", "#8E44AD", "#7F8C8D",
]


def _hex_for_batch(batch: str, registry: dict) -> str:
    if batch not in registry:
        registry[batch] = BATCH_PALETTE[len(registry) % len(BATCH_PALETTE)]
    return registry[batch]


def _safe_float(v) -> float:
    try:
        return float(v) if v not in ("", None, "NA", "nan") else 0.0
    except (TypeError, ValueError):
        return 0.0


def _safe_int(v) -> int:
    try:
        return int(float(v)) if v not in ("", None, "NA", "nan") else 0
    except (TypeError, ValueError):
        return 0


def _read_tsv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _read_csv_any(path: Path) -> list[dict]:
    """Auto-detect tab or comma."""
    if not path.exists():
        return []
    with open(path) as fh:
        head = fh.readline()
        sep = "\t" if "\t" in head else ","
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter=sep))


def _fig_to_b64(fig, dpi: int = 140) -> str:
    """Render a matplotlib Figure to a PNG data: URI."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("ascii")


# ═════════════════════════════════════════════════════════════════════════════
# Data collection
# ═════════════════════════════════════════════════════════════════════════════

@dataclass
class SampleRow:
    """Per-sample joined record used across all sections."""
    sample_id:    str
    platform:     str = ""
    qc_flag:      str = ""
    qc_notes:     str = ""
    sylph_ani:    float = 0.0
    sylph_cov:    float = 0.0
    align_rate:   float = 0.0
    est_depth:    float = 0.0
    contam:       float = 0.0
    raw_illu_reads:  int = 0
    raw_illu_bp:     int = 0
    clean_illu_reads: int = 0
    clean_illu_bp:    int = 0
    clean_nano_bp:    int = 0
    species:      str = ""
    clade:        str = ""
    reference_slug: str = ""
    classification: str = ""
    batch:        str = ""
    is_external:  bool = False
    collection_date: str = ""
    has_amr:      bool = False
    has_cnv:      bool = False

    @property
    def raw_gb(self) -> float:
        return self.raw_illu_bp / 1e9

    @property
    def clean_gb(self) -> float:
        return (self.clean_illu_bp + self.clean_nano_bp) / 1e9

    @property
    def retention_pct(self) -> float:
        return (100.0 * self.clean_illu_bp / self.raw_illu_bp) if self.raw_illu_bp else 0.0

    @property
    def species_label(self) -> str:
        """Concise label for the per-species panel.  C. auris uses the clade."""
        if not self.species:
            return "Unknown"
        if self.species.lower().startswith(("candidozyma auris", "candida auris")):
            if self.clade and len(self.clade) <= 5 and not self.clade.lower().startswith(("not", "n/a", "nan")):
                return f"C. auris clade {self.clade}"
            return "C. auris"
        # Other species: short form "C. albicans"
        parts = self.species.split()
        return f"{parts[0][0]}. {parts[1]}" if len(parts) >= 2 else self.species


def _load_samplesheet(path: Path) -> dict[str, dict]:
    """Indexed by sample id; tolerates missing path."""
    out: dict[str, dict] = {}
    if not path or not path.exists():
        return out
    with open(path) as fh:
        for r in csv.DictReader(fh):
            sid = (r.get("id") or "").strip()
            if sid:
                out[sid] = r
    return out


def _load_fastp_stats(sample_dir: Path) -> tuple[int, int, int, int]:
    """Return (raw_reads, raw_bp, clean_reads, clean_bp) from fastp.json.
    Returns zeros when no JSON is present (e.g. nanopore-only samples).
    """
    cand = list(sample_dir.glob("fastp/*.fastp.json"))
    if not cand:
        return 0, 0, 0, 0
    try:
        with open(cand[0]) as fh:
            data = json.load(fh)
        s = data.get("summary", {})
        b = s.get("before_filtering", {})
        a = s.get("after_filtering",  {})
        return (int(b.get("total_reads", 0)), int(b.get("total_bases", 0)),
                int(a.get("total_reads", 0)), int(a.get("total_bases", 0)))
    except (json.JSONDecodeError, OSError):
        return 0, 0, 0, 0


def collect(results_dir: Path, samplesheet_path: Optional[Path]) -> list[SampleRow]:
    """Walk the results tree and join per-sample data into SampleRow objects."""
    ss = _load_samplesheet(samplesheet_path) if samplesheet_path else {}
    qc_dir = results_dir / "01_qc"
    sp_dir = results_dir / "03_identification"
    amr_dir = results_dir / "06_amr"
    cnv_dir = results_dir / "08_cnv"

    rows: list[SampleRow] = []
    if not qc_dir.exists():
        return rows

    for flag in sorted(qc_dir.glob("*/*.qc_flag.tsv")):
        recs = _read_tsv(flag)
        if not recs:
            continue
        r = recs[0]
        sid = (r.get("sample_id") or flag.parent.name).strip()
        sr = SampleRow(sample_id=sid)
        sr.platform     = (r.get("seq_platform") or "").strip()
        sr.qc_flag      = (r.get("qc_flag") or "").strip()
        sr.qc_notes     = (r.get("qc_notes") or "").strip()
        sr.sylph_ani    = _safe_float(r.get("sylph_ani"))
        sr.sylph_cov    = _safe_float(r.get("sylph_cov_pct"))
        sr.align_rate   = _safe_float(r.get("align_rate_pct"))
        sr.est_depth    = _safe_float(r.get("est_depth_x"))
        sr.contam       = _safe_float(r.get("contam_fraction_pct"))
        sr.clean_illu_bp = _safe_int(r.get("clean_bases_illumina"))
        sr.clean_nano_bp = _safe_int(r.get("clean_bases_nanopore"))

        # Raw + clean read counts pulled from the published fastp JSON. The
        # qc_flag.tsv only carries clean bases so we open the fastp report
        # to populate raw_reads/raw_bp/clean_reads as well.
        rr, rb, cr, cb = _load_fastp_stats(qc_dir / sid)
        sr.raw_illu_reads  = rr
        sr.raw_illu_bp     = rb
        sr.clean_illu_reads = cr
        if cb and not sr.clean_illu_bp:
            sr.clean_illu_bp = cb

        sp_recs = _read_tsv(sp_dir / sid / f"{sid}.species_call.tsv")
        if sp_recs:
            sr.species          = (sp_recs[0].get("species_assigned") or "").strip()
            sr.clade            = (sp_recs[0].get("clade") or "").strip()
            sr.reference_slug   = (sp_recs[0].get("reference_slug") or "").strip()
            sr.classification   = (sp_recs[0].get("classification") or "").strip()

        # Samplesheet enrichment (batch, is_external, collection_date)
        if sid in ss:
            sr.batch           = (ss[sid].get("batch") or "").strip()
            sr.is_external     = (ss[sid].get("is_external") or "").strip().lower() == "true"
            sr.collection_date = (ss[sid].get("collection_date") or "").strip()

        # AMR presence: any non-empty resistance_mutations across the rows
        amr_recs = _read_tsv(amr_dir / sid / f"{sid}.resistance_report.tsv")
        sr.has_amr = any((rr.get("resistance_mutations") or "").strip() for rr in amr_recs)

        # CNV presence: any event != 'normal'
        cnv_recs = _read_tsv(cnv_dir / sid / f"{sid}.cnv_events.tsv")
        sr.has_cnv = any((rr.get("event") or "normal").strip() not in ("normal", "")
                         for rr in cnv_recs)

        rows.append(sr)

    return rows


# ═════════════════════════════════════════════════════════════════════════════
# Branding helpers
# ═════════════════════════════════════════════════════════════════════════════

def load_org_info(branding_dir: Path) -> dict[str, str]:
    info: dict[str, str] = {}
    p = branding_dir / "org.info"
    if not p.exists():
        return info
    for line in p.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, _, v = line.partition("=")
        info[k.strip()] = v.strip()
    return info


def encode_logo_b64(branding_dir: Path) -> str:
    p = branding_dir / "logo.png"
    if not p.exists():
        return ""
    return "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode("ascii")


def file_to_b64(path: Path, mime: str = "image/png") -> str:
    if not path.exists():
        return ""
    return f"data:{mime};base64," + base64.b64encode(path.read_bytes()).decode("ascii")


# ═════════════════════════════════════════════════════════════════════════════
# CSS
# ═════════════════════════════════════════════════════════════════════════════

CSS = f"""
*{{box-sizing:border-box;margin:0;padding:0;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;}}
body{{background:{BG};color:{DARK};line-height:1.45;}}
.wrap{{max-width:1320px;margin:0 auto;padding:24px;}}

.banner{{
  display:flex;align-items:center;gap:24px;
  background:linear-gradient(135deg,{BLUE} 0%,{LBLUE} 100%);
  color:#fff;padding:24px 28px;border-radius:14px;
  box-shadow:0 6px 24px rgba(44,95,138,0.18);
  margin-bottom:24px;
}}
.banner img{{background:#fff;padding:8px 12px;border-radius:8px;}}
.banner img.logo-epicandi{{height:78px;}}
.banner img.logo-org     {{height:56px;margin-left:-4px;}}
.banner .titles{{flex:1;}}
.banner h1{{font-size:22px;font-weight:700;margin-bottom:4px;}}
.banner h2{{font-size:14px;font-weight:400;opacity:0.92;}}
.banner .org{{font-size:12px;opacity:0.88;text-align:right;line-height:1.6;}}
.banner .org a{{color:#fff;text-decoration:underline;}}

.section{{background:#fff;border-radius:12px;padding:20px 24px;box-shadow:0 2px 8px rgba(0,0,0,0.04);margin-bottom:20px;}}
.section h3{{font-size:16px;font-weight:700;color:{BLUE};margin-bottom:14px;border-left:4px solid {ORANGE};padding-left:10px;}}
.section h4{{font-size:13px;font-weight:600;color:{DARK};margin:14px 0 8px;}}

.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:16px;}}
.card{{background:#fff;border-radius:10px;padding:18px;box-shadow:0 1px 4px rgba(0,0,0,0.06);border-top:3px solid {BLUE};}}
.card.platform{{border-top-color:{LBLUE};}}
.card.pass    {{border-top-color:{GREEN};}}
.card.fail    {{border-top-color:{RED};}}
.card.amr     {{border-top-color:{ORANGE};}}
.card.cnv     {{border-top-color:#9B59B6;}}
.card .label{{font-size:11px;text-transform:uppercase;letter-spacing:0.04em;color:{GREY};margin-bottom:6px;}}
.card .big  {{font-size:30px;font-weight:700;color:{DARK};line-height:1.1;}}
.card .sub  {{font-size:12px;color:{GREY};margin-top:4px;}}

.platform-grid{{display:grid;grid-template-columns:repeat(3,1fr);gap:8px;margin-top:10px;}}
.platform-grid div{{background:{BG};padding:8px;border-radius:6px;text-align:center;}}
.platform-grid .pct{{font-size:18px;font-weight:700;color:{BLUE};}}
.platform-grid .tag{{font-size:10px;color:{GREY};text-transform:uppercase;letter-spacing:0.04em;}}

img.fig{{max-width:100%;height:auto;display:block;margin:8px auto;border-radius:6px;}}

table.dat{{border-collapse:collapse;width:100%;font-size:12px;margin-top:8px;}}
table.dat th{{background:{BG};text-align:left;padding:6px 10px;color:{DARK};font-weight:600;border-bottom:2px solid {BLUE};}}
table.dat td{{padding:5px 10px;border-bottom:1px solid #ECEEF0;}}
table.dat tr:hover{{background:#FAFBFC;}}
.pill{{display:inline-block;padding:2px 8px;border-radius:10px;font-size:11px;font-weight:600;color:#fff;}}
.pill.pass{{background:{GREEN};}}
.pill.fail{{background:{RED};}}
.pill.warn{{background:{ORANGE};}}

/* Heat-bar fill inside numeric cells (proportional to value within its column). */
td.barcell{{position:relative;}}
td.barcell .fill{{position:absolute;top:2px;bottom:2px;left:4px;border-radius:3px;opacity:0.32;z-index:0;}}
td.barcell .val{{position:relative;z-index:1;font-variant-numeric:tabular-nums;}}

.batch-strip{{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:6px;vertical-align:middle;}}
.legend{{display:flex;flex-wrap:wrap;gap:14px;margin-top:8px;font-size:11px;color:{GREY};}}
.legend span{{display:inline-flex;align-items:center;}}
.legend .sw{{display:inline-block;width:12px;height:12px;border-radius:3px;margin-right:5px;}}

.footer{{font-size:11px;color:{GREY};text-align:center;padding:18px 0;}}
.muted{{color:{GREY};font-size:12px;}}
.two-col{{display:grid;grid-template-columns:2fr 1fr;gap:20px;align-items:start;}}
"""


# ═════════════════════════════════════════════════════════════════════════════
# Section renderers
# ═════════════════════════════════════════════════════════════════════════════

def render_banner(run_name: str, org: dict[str, str],
                  logo_uri: str, epicandi_uri: str = "") -> str:
    org_html = ""
    if org:
        org_html = f"""
            <div class="org">
              <strong>{org.get('name','')}</strong><br>
              {org.get('department','')}<br>
              {org.get('address','')}<br>
              <a href="mailto:{org.get('contact','')}">{org.get('contact','')}</a>
              &middot; <a href="{org.get('web','')}">{org.get('web','').replace('https://','').replace('http://','')}</a>
            </div>"""
    # Epicandi pipeline logo first (left, larger), then the EPIMOL org logo.
    epicandi_html = (f'<img src="{epicandi_uri}" alt="epicandi" class="logo-epicandi">'
                     if epicandi_uri else "")
    logo_html     = (f'<img src="{logo_uri}" alt="logo" class="logo-org">'
                     if logo_uri else "")
    return f"""
    <header class="banner">
      {epicandi_html}{logo_html}
      <div class="titles">
        <h1>epicandi — Run report</h1>
        <h2>{run_name}  &middot;  Generated {datetime.now().strftime('%Y-%m-%d %H:%M')}</h2>
      </div>
      {org_html}
    </header>
    """


def render_summary_cards(rows: list[SampleRow]) -> str:
    total = len(rows)
    short = sum(1 for r in rows if r.platform == "illumina")
    long_ = sum(1 for r in rows if r.platform == "nanopore")
    hybrid = sum(1 for r in rows if r.platform == "hybrid")
    pass_n = sum(1 for r in rows if r.qc_flag == "PASS")
    fail_n = total - pass_n
    amr_n = sum(1 for r in rows if r.has_amr)
    cnv_n = sum(1 for r in rows if r.has_cnv)
    clean_gb = sum(r.clean_gb for r in rows)

    pct = lambda n: (n / total * 100.0) if total else 0.0

    platform_card = f"""
    <div class="card platform">
      <div class="label">Samples analysed</div>
      <div class="big">{total}</div>
      <div class="platform-grid">
        <div><div class="pct">{short}</div><div class="tag">short</div></div>
        <div><div class="pct">{long_}</div><div class="tag">long</div></div>
        <div><div class="pct">{hybrid}</div><div class="tag">hybrid</div></div>
      </div>
    </div>"""

    gb_card = f"""
    <div class="card">
      <div class="label">Sequencing volume</div>
      <div class="big">{clean_gb:.1f}<span style="font-size:14px;color:{GREY};font-weight:400;"> GB clean</span></div>
      <div class="sub">raw not joined yet &middot; clean used as proxy</div>
    </div>"""

    pass_card = f"""
    <div class="card pass">
      <div class="label">QC pass</div>
      <div class="big">{pass_n}<span style="font-size:14px;color:{GREY};font-weight:400;"> / {total}</span></div>
      <div class="sub">{pct(pass_n):.1f}%  &middot;  {fail_n} failed</div>
    </div>"""

    amr_card = f"""
    <div class="card amr">
      <div class="label">Samples with AMR mutations</div>
      <div class="big">{amr_n}</div>
      <div class="sub">{pct(amr_n):.1f}% of cohort</div>
    </div>"""

    cnv_card = f"""
    <div class="card cnv">
      <div class="label">Samples with CNV events</div>
      <div class="big">{cnv_n}</div>
      <div class="sub">{pct(cnv_n):.1f}% of cohort</div>
    </div>"""

    return f"""
    <section class="section">
      <h3>Run summary</h3>
      <div class="cards">
        {platform_card}{gb_card}{pass_card}{amr_card}{cnv_card}
      </div>
    </section>
    """


# ─── Phase 2 — species composition + identification chart (Plotly) ─────────
#
# These two charts are rendered with Plotly.js (loaded once from CDN in the
# document head) so the user can hover, toggle batches and export PNG/SVG from
# the browser.  Each bar carries a numeric annotation so that values stay
# readable when the report is printed to PDF.

_PLOTLY_DIV_COUNTER = 0


def _plotly_div(traces: list[dict], layout: dict, height_px: int = 380) -> str:
    """Return a self-contained <div>+<script> pair that draws a Plotly chart.
    The plot call is wrapped in try/catch so that a failure in one chart
    cannot block the others or break the rest of the page.

    Compact charts (height ≤ 200 px — the mini-histograms inside the summary
    cards) get no `min-width`, so they shrink to their card width.  Larger
    charts get a `min-width` to guard against a still-zero-width container at
    first render in some browsers.
    """
    global _PLOTLY_DIV_COUNTER
    _PLOTLY_DIV_COUNTER += 1
    div_id = f"plotly_div_{_PLOTLY_DIV_COUNTER}"
    is_compact = height_px <= 200
    min_w = "" if is_compact else "min-width:600px;"
    return (
        f'<div id="{div_id}" style="width:100%;{min_w}height:{height_px}px;'
        f'box-sizing:content-box;"></div>'
        f'<script>try{{Plotly.newPlot("{div_id}",'
        f' {json.dumps(traces)}, {json.dumps(layout)},'
        f' {{responsive:true, displaylogo:false, modeBarButtonsToRemove:["lasso2d","select2d"]}})'
        f'  .then(function(){{ '
        f'    /* Force a final relayout once the page has settled so charts '
        f'       whose container width was 0 at first render get re-painted. */'
        f'    requestAnimationFrame(function(){{ try{{Plotly.Plots.resize("{div_id}");}}catch(_){{}}; }});'
        f'  }});'
        f'}}catch(_e){{document.getElementById("{div_id}").innerHTML='
        f'"<div style=\\"padding:14px;color:#c00;border:1px dashed #c00;\\">'
        f'Plotly render error: "+_e.message+"</div>";}}</script>'
    )


def _bar_total_labels(totals: list[int]) -> list[str]:
    """Render counts shown above each bar (empty string for zeros)."""
    return [str(t) if t > 0 else "" for t in totals]


def _bar_cell(val_str: str, val: float, vmin: float, vmax: float,
              color_hex: str) -> str:
    """Return a <td class="barcell"> whose background bar grows with the value.
    Used for numeric columns in the per-sample QC / assembly tables so the
    reader can compare rows at a glance.
    """
    if vmax <= vmin or val <= 0:
        pct = 0.0
    else:
        pct = max(0.0, min(100.0, 100.0 * (val - vmin) / (vmax - vmin)))
    return (f'<td class="barcell">'
            f'<span class="fill" style="background:{color_hex};width:{pct:.1f}%;"></span>'
            f'<span class="val">{val_str}</span></td>')


# Column → colour group used for cell-bar shading.
COLOUR_FASTQ    = BLUE        # raw_reads / clean_reads / GB / retention
COLOUR_SYLPH    = ORANGE      # sylph ANI / cov / contam / align / depth
COLOUR_ASSEMBLY = GREEN       # QUAST + BUSCO numeric metrics


def render_species_batch_section(rows: list[SampleRow],
                                 batch_colors: dict[str, str]) -> str:
    """Interactive stacked bar per species detected (Plotly).

    X-axis: collection year if any sample has collection_date, otherwise batch.
    Number of samples is annotated on top of each bar.
    Excludes is_external=true.
    """
    eligible = [r for r in rows if r.qc_flag == "PASS" and not r.is_external]
    if not eligible:
        return ""

    by_species: dict[str, list[SampleRow]] = defaultdict(list)
    for r in eligible:
        by_species[r.species_label].append(r)

    has_dates = any(r.collection_date for r in eligible)
    batches_all = sorted({r.batch for r in eligible if r.batch})
    show_legend = len(batches_all) > 1                    # one or zero batches: hide legend.

    # Shared X range + bar width across species panels.  Without this, Plotly
    # auto-fits each subplot to its own subset of years — so tropicalis or
    # parapsilosis with only one collection year render a single bar that
    # spans the entire plot (the "huge mountain" effect).  Locking the range
    # to the whole-cohort span and giving every bar an explicit width of 0.6
    # years makes the visual weight of a "1 sample / 1 year" species
    # comparable to a "60 sample / 4 year" species.
    #
    # The Y axis is also locked to the same range across every species panel
    # so a bar of height 1 (e.g. tropicalis, one isolate) doesn't visually
    # match a bar of height 9 (auris in a busy year) just because Plotly
    # auto-scales each panel independently.
    BAR_WIDTH_YEARS = 0.6
    xaxis_range = None
    yaxis_max = None
    if has_dates:
        years_all_global = sorted({int(r.collection_date[:4])
                                   for r in eligible if r.collection_date})
        if years_all_global:
            ymin, ymax = years_all_global[0], years_all_global[-1]
            if ymin == ymax:
                ymin -= 1; ymax += 1
            xaxis_range = [ymin - 0.5, ymax + 0.5]

            # Compute the global max stacked total per (species, year).  The
            # annotation labels sit yshift:8 above each bar, so we add
            # headroom (≥1, ≥15 % of max) so they don't get clipped.
            max_total = 0
            for _sp, _subset in by_species.items():
                for y in {int(r.collection_date[:4])
                          for r in _subset if r.collection_date}:
                    tot = sum(1 for r in _subset
                              if r.collection_date and int(r.collection_date[:4]) == y)
                    if tot > max_total:
                        max_total = tot
            if max_total:
                yaxis_max = max_total + max(1, int(round(max_total * 0.15)))

    chart_blocks: list[str] = []
    for sp_label, subset in sorted(by_species.items(), key=lambda kv: -len(kv[1])):

        if has_dates and batches_all:
            # X = collection year as integer, stacked by batch.
            # Using `type: "linear"` with int x (not category strings) — works
            # around a Plotly bug where 4-digit year strings under category mode
            # render the axes + annotations but skip drawing the bars.
            years = sorted({int(r.collection_date[:4]) for r in subset if r.collection_date})
            traces = []
            totals = [0] * len(years)
            for b in batches_all:
                ys = [sum(1 for r in subset
                          if r.batch == b and int(r.collection_date[:4]) == y) for y in years]
                if sum(ys) == 0:
                    continue
                totals = [t + v for t, v in zip(totals, ys)]
                traces.append({
                    "type": "bar", "name": b, "x": years, "y": ys,
                    "width": [BAR_WIDTH_YEARS] * len(years),
                    "marker": {"color": _hex_for_batch(b, batch_colors),
                               "line": {"color": "white", "width": 0.5}},
                    "hovertemplate": f"<b>{b}</b><br>%{{x}}: %{{y}}<extra></extra>",
                })
            annotations = [{
                "x": yr, "y": tot, "text": str(tot),
                "showarrow": False, "yshift": 8,
                "font": {"size": 11, "color": DARK},
            } for yr, tot in zip(years, totals) if tot > 0]
            xaxis_cfg = {"title": "Collection year", "type": "linear",
                         "dtick": 1, "tickformat": "d"}
            if xaxis_range is not None:
                xaxis_cfg["range"] = xaxis_range
            yaxis_cfg = {"title": "Samples"}
            if yaxis_max is not None:
                yaxis_cfg["range"] = [0, yaxis_max]
                yaxis_cfg["dtick"] = 1 if yaxis_max <= 12 else None
            layout = {
                "title": {"text": f"{sp_label} — {len(subset)} samples",
                          "font": {"size": 14, "color": BLUE}},
                "barmode": "stack",
                "xaxis": xaxis_cfg,
                "yaxis": yaxis_cfg,
                "showlegend": show_legend,
                "legend": {"orientation": "v", "x": 1.02, "y": 1.0,
                           "xanchor": "left", "yanchor": "top",
                           "title": {"text": "Batch"}},
                "margin": {"l": 60, "r": 160 if show_legend else 40,
                           "t": 50, "b": 50},
                "annotations": annotations,
                "plot_bgcolor": "white", "paper_bgcolor": "white",
            }
            chart_blocks.append(_plotly_div(traces, layout, height_px=420))

        elif batches_all:
            # No dates: one bar per batch.
            counts = [sum(1 for r in subset if r.batch == b) for b in batches_all]
            traces = [{
                "type": "bar", "x": batches_all, "y": counts,
                "marker": {"color": [_hex_for_batch(b, batch_colors)
                                     for b in batches_all]},
                "text": _bar_total_labels(counts), "textposition": "outside",
                "hovertemplate": "%{x}: %{y}<extra></extra>",
                "showlegend": False,
            }]
            layout = {
                "title": {"text": f"{sp_label} — {len(subset)} samples",
                          "font": {"size": 14, "color": BLUE}},
                "xaxis": {"title": "Batch", "type": "category"},
                "yaxis": {"title": "Samples"},
                "showlegend": False,
                "margin": {"l": 60, "r": 40, "t": 50, "b": 60},
                "plot_bgcolor": "white", "paper_bgcolor": "white",
            }
            chart_blocks.append(_plotly_div(traces, layout, height_px=380))
        else:
            # No batch metadata at all: total only.
            traces = [{
                "type": "bar", "x": [sp_label], "y": [len(subset)],
                "marker": {"color": BLUE},
                "text": [str(len(subset))], "textposition": "outside",
                "showlegend": False,
            }]
            layout = {
                "title": {"text": f"{sp_label} — {len(subset)} samples",
                          "font": {"size": 14, "color": BLUE}},
                "yaxis": {"title": "Samples"},
                "showlegend": False, "xaxis": {"showticklabels": False},
                "margin": {"l": 60, "r": 40, "t": 50, "b": 30},
                "plot_bgcolor": "white", "paper_bgcolor": "white",
            }
            chart_blocks.append(_plotly_div(traces, layout, height_px=280))

    if not chart_blocks:
        return ""

    note_bits = ["Excluding samples marked as external."]
    if has_dates and batches_all:
        note_bits.append("Stacked by batch across collection years.")
    elif batches_all:
        note_bits.append("Batch breakdown (no collection_date in samplesheet).")
    note = f'<p class="muted">{" ".join(note_bits)}</p>'

    return f"""
    <section class="section">
      <h3>Sample composition per species</h3>
      {note}
      {''.join(chart_blocks)}
    </section>
    """


def render_identification_chart(rows: list[SampleRow],
                                batch_colors: dict[str, str]) -> str:
    """Interactive stacked bar of PASS samples per reference (Plotly).

    - References sorted by descending total count.
    - When more than one batch is present, stacks are coloured by batch and the
      legend is moved outside the plot area.  With a single batch (or none) the
      bars share the corporate blue and no legend is shown.
    - Each bar carries the total count above it for print/PDF readability.
    Excludes is_external=true.
    """
    eligible = [r for r in rows if r.qc_flag == "PASS" and not r.is_external]
    if not eligible:
        return ""

    refs = sorted({r.reference_slug or "unknown" for r in eligible},
                  key=lambda x: -sum(1 for rr in eligible
                                     if (rr.reference_slug or "unknown") == x))
    batches_all = sorted({r.batch for r in eligible if r.batch})
    multi_batch = len(batches_all) > 1

    totals = [sum(1 for r in eligible if (r.reference_slug or "unknown") == ref)
              for ref in refs]

    if multi_batch:
        traces = []
        for b in batches_all:
            ys = [sum(1 for r in eligible
                      if r.batch == b and (r.reference_slug or "unknown") == ref)
                  for ref in refs]
            if sum(ys) == 0:
                continue
            traces.append({
                "type": "bar", "name": b, "x": refs, "y": ys,
                "marker": {"color": _hex_for_batch(b, batch_colors),
                           "line": {"color": "white", "width": 0.5}},
                "hovertemplate": f"<b>{b}</b><br>%{{x}}: %{{y}}<extra></extra>",
            })
        annotations = [{
            "x": ref, "y": tot, "text": str(tot),
            "showarrow": False, "yshift": 10,
            "font": {"size": 11, "color": DARK},
        } for ref, tot in zip(refs, totals) if tot > 0]
        layout = {
            "title": {"text": "Identification — reference hits",
                      "font": {"size": 14, "color": BLUE}},
            "barmode": "stack",
            "xaxis": {"title": "Reference", "type": "category", "tickangle": -30},
            "yaxis": {"title": "Samples (PASS)"},
            "showlegend": True,
            "legend": {"orientation": "v", "x": 1.02, "y": 1.0,
                       "xanchor": "left", "yanchor": "top",
                       "title": {"text": "Batch"}},
            "margin": {"l": 60, "r": 180, "t": 50, "b": 110},
            "annotations": annotations,
            "plot_bgcolor": "white", "paper_bgcolor": "white",
        }
    else:
        traces = [{
            "type": "bar", "x": refs, "y": totals,
            "marker": {"color": (_hex_for_batch(batches_all[0], batch_colors)
                                 if batches_all else BLUE)},
            "text": _bar_total_labels(totals), "textposition": "outside",
            "hovertemplate": "%{x}: %{y}<extra></extra>",
            "showlegend": False,
        }]
        layout = {
            "title": {"text": "Identification — reference hits",
                      "font": {"size": 14, "color": BLUE}},
            "xaxis": {"title": "Reference", "type": "category", "tickangle": -30},
            "yaxis": {"title": "Samples (PASS)"},
            "showlegend": False,
            "margin": {"l": 60, "r": 40, "t": 50, "b": 110},
            "plot_bgcolor": "white", "paper_bgcolor": "white",
        }

    chart_html = _plotly_div(traces, layout,
                             height_px=max(380, 60 + 18 * len(refs) + 100))

    # Mini-table beside the chart.
    table_rows = "".join(f"<tr><td>{ref}</td><td>{tot}</td></tr>"
                         for ref, tot in zip(refs, totals))
    table = ('<table class="dat"><thead><tr><th>Reference</th>'
             '<th>Identified (PASS)</th></tr></thead><tbody>'
             f'{table_rows}</tbody></table>')

    return f"""
    <section class="section">
      <h3>Identification — per reference</h3>
      <p class="muted">Excluding samples marked as external. {len(eligible)} PASS samples shown.</p>
      <div class="two-col">
        <div>{chart_html}</div>
        <div>{table}</div>
      </div>
    </section>
    """


# ─── Phase 3 — Fastq QC ────────────────────────────────────────────────────

def render_fastq_qc(rows: list[SampleRow],
                    batch_colors: dict[str, str],
                    secondary_html_path: Optional[Path]) -> str:
    if not rows:
        return ""

    # ── Run-level totals: 5 cards, the Reads / GBs / Retention / Depth cards
    # carry a small Plotly histogram inside so the user can see the per-sample
    # distribution and hover over individual bars.
    total_raw_bp   = sum(r.raw_illu_bp   for r in rows)
    total_clean_bp = sum(r.clean_illu_bp + r.clean_nano_bp for r in rows)
    total_raw_reads   = sum(r.raw_illu_reads   for r in rows)
    total_clean_reads = sum(r.clean_illu_reads for r in rows)
    retention = (100.0 * total_clean_bp / total_raw_bp) if total_raw_bp else 0.0
    pass_n = sum(1 for r in rows if r.qc_flag == "PASS")
    median_retention = float(np.median([r.retention_pct for r in rows
                                        if r.retention_pct > 0])) if rows else 0.0
    median_depth     = float(np.median([r.est_depth for r in rows
                                        if r.est_depth > 0])) if rows else 0.0

    def _mini_hist(traces: list[dict], shapes: list[dict] = None,
                   showlegend: bool = False, height: int = 110) -> str:
        layout = {
            "margin": {"l": 30, "r": 8, "t": 6, "b": 22},
            "xaxis": {"tickfont": {"size": 9}, "showgrid": False},
            "yaxis": {"tickfont": {"size": 9}, "showgrid": False, "title": ""},
            "showlegend": showlegend,
            "legend": {"orientation": "h", "x": 0, "y": 1.15,
                       "font": {"size": 9}, "yanchor": "bottom"},
            "barmode": "overlay",
            "plot_bgcolor": "white", "paper_bgcolor": "white",
            "shapes": shapes or [],
        }
        return _plotly_div(traces, layout, height_px=height)

    # Reads card — overlaid RAW (red) + CLEAN (blue) histograms.
    reads_traces = [
        {"type": "histogram", "x": [r.raw_illu_reads/1e6 for r in rows],
         "name": "raw", "marker": {"color": RED},   "opacity": 0.55,
         "hovertemplate": "raw: %{x:.1f} M reads<br>n=%{y}<extra></extra>"},
        {"type": "histogram", "x": [r.clean_illu_reads/1e6 for r in rows],
         "name": "clean", "marker": {"color": BLUE}, "opacity": 0.55,
         "hovertemplate": "clean: %{x:.1f} M reads<br>n=%{y}<extra></extra>"},
    ]
    reads_card = (
        '<div class="card platform" style="grid-column:span 2;">'
        '<div class="label">Reads (M, distribution per sample)</div>'
        f'<div class="big" style="font-size:18px;">'
        f'{total_raw_reads/1e6:.1f} M <span style="font-size:11px;color:{GREY};">raw</span> · '
        f'{total_clean_reads/1e6:.1f} M <span style="font-size:11px;color:{GREY};">clean</span></div>'
        + _mini_hist(reads_traces, showlegend=True, height=130)
        + '</div>')

    # GBs card — overlaid RAW (red) + CLEAN (blue) GB histograms.
    gb_traces = [
        {"type": "histogram", "x": [r.raw_gb   for r in rows],
         "name": "raw",  "marker": {"color": RED},  "opacity": 0.55,
         "hovertemplate": "raw: %{x:.2f} GB<br>n=%{y}<extra></extra>"},
        {"type": "histogram", "x": [r.clean_gb for r in rows],
         "name": "clean","marker": {"color": BLUE}, "opacity": 0.55,
         "hovertemplate": "clean: %{x:.2f} GB<br>n=%{y}<extra></extra>"},
    ]
    gb_card = (
        '<div class="card platform" style="grid-column:span 2;">'
        '<div class="label">Gigabases (distribution per sample)</div>'
        f'<div class="big" style="font-size:18px;">'
        f'{total_raw_bp/1e9:.1f} GB <span style="font-size:11px;color:{GREY};">raw</span> · '
        f'{total_clean_bp/1e9:.1f} GB <span style="font-size:11px;color:{GREY};">clean</span></div>'
        + _mini_hist(gb_traces, showlegend=True, height=130)
        + '</div>')

    # Retention card — single histogram.
    ret_traces = [{
        "type": "histogram", "x": [r.retention_pct for r in rows if r.retention_pct > 0],
        "marker": {"color": GREEN}, "opacity": 0.75,
        "hovertemplate": "%{x:.1f} % · n=%{y}<extra></extra>",
    }]
    ret_card = (
        '<div class="card pass">'
        '<div class="label">% Retention (per sample)</div>'
        f'<div class="big" style="font-size:18px;">{retention:.1f} %</div>'
        f'<div class="muted" style="font-size:10px;">median {median_retention:.1f} %</div>'
        + _mini_hist(ret_traces, height=110) + '</div>')

    # Helper for the single-histogram identification cards (ANI, cov, align,
    # depth) — orange family with a dashed threshold line.
    def _ident_card(label_top, value_top, sub, vals, threshold,
                    hover_unit, value_fmt="{:.1f}"):
        shapes = []
        if threshold is not None:
            shapes = [{
                "type": "line", "xref": "x", "yref": "paper",
                "x0": threshold, "x1": threshold, "y0": 0, "y1": 1,
                "line": {"color": GREY, "width": 1, "dash": "dash"},
            }]
        traces = [{
            "type": "histogram", "x": vals,
            "marker": {"color": ORANGE}, "opacity": 0.75,
            "hovertemplate": f"%{{x:{'.1f' if value_fmt!='{:.0f}' else '.0f'}}} {hover_unit}<br>n=%{{y}}<extra></extra>",
        }]
        return (
            '<div class="card amr">'
            f'<div class="label">{label_top}</div>'
            f'<div class="big" style="font-size:18px;">{value_top}</div>'
            f'<div class="muted" style="font-size:10px;">{sub}</div>'
            + _mini_hist(traces, shapes, height=110) + '</div>')

    median_ani   = float(np.median([r.sylph_ani  for r in rows if r.sylph_ani  > 0])) if rows else 0.0
    median_cov   = float(np.median([r.sylph_cov  for r in rows if r.sylph_cov  > 0])) if rows else 0.0
    median_align = float(np.median([r.align_rate for r in rows if r.align_rate > 0])) if rows else 0.0

    ani_card   = _ident_card("Sylph ANI (%)", f"{median_ani:.2f} %",
                             "median &middot; dashed = 97% threshold",
                             [r.sylph_ani for r in rows if r.sylph_ani > 0],
                             97.0, "%")
    cov_card   = _ident_card("Sylph coverage (%)", f"{median_cov:.1f} %",
                             "median &middot; dashed = 85% threshold",
                             [r.sylph_cov for r in rows if r.sylph_cov > 0],
                             85.0, "%")
    align_card = _ident_card("Read alignment rate (%)", f"{median_align:.1f} %",
                             "median &middot; dashed = 90% threshold",
                             [r.align_rate for r in rows if r.align_rate > 0],
                             90.0, "%")

    # Mean depth estimated card — single histogram with threshold line at 30x.
    depth_vals = [r.est_depth for r in rows if r.est_depth > 0]
    depth_traces = [{
        "type": "histogram", "x": depth_vals,
        "marker": {"color": ORANGE}, "opacity": 0.75,
        "hovertemplate": "%{x:.0f}×<br>n=%{y}<extra></extra>",
    }]
    depth_shapes = [{
        "type": "line", "xref": "x", "yref": "paper",
        "x0": 30, "x1": 30, "y0": 0, "y1": 1,
        "line": {"color": GREY, "width": 1, "dash": "dash"},
    }]
    depth_card = (
        '<div class="card amr">'
        '<div class="label">Mean depth estimated (×)</div>'
        f'<div class="big" style="font-size:18px;">{median_depth:.0f}×</div>'
        f'<div class="muted" style="font-size:10px;">median &middot; dashed = 30× threshold</div>'
        + _mini_hist(depth_traces, depth_shapes, height=110) + '</div>')

    # Samples PASS card — text only.
    pass_card = (
        '<div class="card pass">'
        '<div class="label">Samples (PASS)</div>'
        f'<div class="big" style="font-size:24px;">{pass_n} <span style="font-size:14px;color:{GREY};">/ {len(rows)}</span></div>'
        f'<div class="muted" style="font-size:11px;">{100.0*pass_n/len(rows):.1f} %</div>'
        '</div>')

    # First row: Samples PASS + Reads (2 cols) + GBs (2 cols) + Retention.
    # Second row: ANI / cov / align / depth (4 single-column cards).
    row1 = (f'<div class="cards" style="grid-template-columns:repeat(6,1fr);">'
            f'{pass_card}{reads_card}{gb_card}{ret_card}</div>')
    row2 = (f'<div class="cards" style="grid-template-columns:repeat(4,1fr);margin-top:14px;">'
            f'{ani_card}{cov_card}{align_card}{depth_card}</div>')
    stats_html = row1 + row2

    # ── Plotly 4-panel: ANI / cov / align / depth with hover on each point.
    THRESHOLDS = {"sylph_ani": 97.0, "sylph_cov": 85.0,
                  "align_rate": 90.0, "est_depth": 30.0}
    PANEL_TITLE = {
        "sylph_ani":  "Sylph ANI (%)",
        "sylph_cov":  "Sylph coverage (%)",
        "align_rate": "Read alignment rate (%)",
        "est_depth":  "Estimated mean depth (×)",
    }
    PANEL_AXIS = {"sylph_ani": ("x", "y"),  "sylph_cov": ("x2", "y2"),
                  "align_rate": ("x3", "y3"), "est_depth": ("x4", "y4")}

    sids_all = [r.sample_id for r in rows]
    traces, annotations, shapes = [], [], []
    for key in ["sylph_ani", "sylph_cov", "align_rate", "est_depth"]:
        vals = [getattr(r, key) for r in rows]
        cols = [GREEN if v >= THRESHOLDS[key] else RED for v in vals]
        xaxis, yaxis = PANEL_AXIS[key]
        traces.append({
            "type": "scatter", "mode": "markers",
            "x": list(range(len(vals))), "y": vals,
            "text": sids_all,
            "hovertemplate": f"<b>%{{text}}</b><br>{PANEL_TITLE[key]}: %{{y:.2f}}<br>threshold: {THRESHOLDS[key]:g}<extra></extra>",
            "marker": {"color": cols, "size": 7, "line": {"width": 0}},
            "xaxis": xaxis, "yaxis": yaxis, "showlegend": False,
        })
        shapes.append({
            "type": "line", "xref": xaxis, "yref": yaxis,
            "x0": 0, "x1": max(0, len(vals) - 1),
            "y0": THRESHOLDS[key], "y1": THRESHOLDS[key],
            "line": {"color": GREY, "width": 1, "dash": "dash"},
        })
        annotations.append({
            "xref": f"{xaxis} domain", "yref": yaxis,
            "x": 0.99, "y": THRESHOLDS[key], "xanchor": "right", "yanchor": "bottom",
            "text": f" threshold = {THRESHOLDS[key]:g}",
            "showarrow": False, "font": {"size": 10, "color": GREY},
        })
    layout_panel = {
        "grid": {"rows": 2, "columns": 2, "pattern": "independent"},
        "xaxis":  {"showticklabels": False, "title": ""},
        "xaxis2": {"showticklabels": False, "title": ""},
        "xaxis3": {"showticklabels": False, "title": ""},
        "xaxis4": {"showticklabels": False, "title": ""},
        "yaxis":  {"title": PANEL_TITLE["sylph_ani"]},
        "yaxis2": {"title": PANEL_TITLE["sylph_cov"]},
        "yaxis3": {"title": PANEL_TITLE["align_rate"]},
        "yaxis4": {"title": PANEL_TITLE["est_depth"]},
        "shapes": shapes, "annotations": annotations,
        "margin": {"l": 60, "r": 30, "t": 50, "b": 30},
        "plot_bgcolor": "white", "paper_bgcolor": "white",
        "title": {"text": "Identification & coverage QC",
                  "font": {"size": 14, "color": BLUE}},
        "showlegend": False,
    }
    panel_html = _plotly_div(traces, layout_panel, height_px=620)

    # ── FAIL minitable (unchanged).
    fails = [r for r in rows if r.qc_flag != "PASS"]
    fail_rows = ""
    for r in fails:
        cells = "".join(f"<td>{v}</td>" for v in [
            r.sample_id, r.reference_slug or r.species, f"{r.sylph_cov:.1f}",
            f"{r.align_rate:.1f}", f"{r.est_depth:.1f}"
        ])
        pill_class = "pass" if r.qc_flag == "PASS" else "fail"
        fail_rows += (f"<tr>{cells}"
                      f"<td><span class='pill {pill_class}'>{r.qc_flag}</span></td>"
                      f"<td class='muted'>{r.qc_notes}</td></tr>")
    fail_table = ("<table class='dat'><thead><tr>"
                  "<th>Sample</th><th>match_sylph</th><th>%cov</th>"
                  "<th>%align</th><th>depth</th><th>Status</th><th>Reason</th>"
                  f"</tr></thead><tbody>{fail_rows or '<tr><td colspan=7 class=muted>(no failures)</td></tr>'}</tbody></table>")

    link_html = ""
    if secondary_html_path:
        rel = secondary_html_path.name
        link_html = (f'<p><a href="{rel}" target="_blank" '
                     f'style="color:{BLUE};font-weight:600;">Open full per-sample QC table →</a></p>')

    return f"""
    <section class="section">
      <h3>Fastq QC</h3>
      {link_html}
      <h4>Run totals</h4>
      {stats_html}
      <h4>Identification &amp; coverage QC</h4>
      {panel_html}
      <h4>Samples failing QC</h4>
      {fail_table}
    </section>
    """


SORT_JS = r"""
<script>
(function(){
  // Numeric value extraction: digit tolerant of commas, %, ×, " M", etc.
  function cellNum(td){
    var raw = (td.querySelector('.val') ? td.querySelector('.val').textContent : td.textContent) || '';
    raw = raw.trim().replace(/,/g,'').replace(/[%×*xX]/g,'').replace(/\s+M$/,'e6').replace(/\s+/g,'');
    var n = parseFloat(raw);
    return isNaN(n) ? null : n;
  }
  function attachSort(table){
    var ths = table.tHead.rows[0].cells;
    Array.from(ths).forEach(function(th, idx){
      th.style.cursor = 'pointer';
      th.title = 'Click to sort';
      var dir = 1;  // 1 ascending, -1 descending
      th.addEventListener('click', function(){
        Array.from(ths).forEach(function(o){o.dataset.sort='';});
        th.dataset.sort = (dir===1?'asc':'desc');
        var tbody = table.tBodies[0];
        var rows  = Array.from(tbody.rows);
        rows.sort(function(a,b){
          var va = cellNum(a.cells[idx]), vb = cellNum(b.cells[idx]);
          if (va===null && vb===null)
            return a.cells[idx].textContent.localeCompare(b.cells[idx].textContent) * dir;
          if (va===null) return 1;
          if (vb===null) return -1;
          return (va-vb) * dir;
        });
        rows.forEach(function(r){tbody.appendChild(r);});
        dir = -dir;
      });
    });
  }
  document.querySelectorAll('table.sortable').forEach(attachSort);
})();
</script>
<style>
table.sortable th{position:relative;}
table.sortable th[data-sort=asc]::after  {content:' ▲';font-size:9px;color:#888;}
table.sortable th[data-sort=desc]::after {content:' ▼';font-size:9px;color:#888;}
</style>
"""

FILTER_JS_TEMPLATE = r"""
<style>
.filter-panel{background:#fff;border-radius:10px;padding:14px 18px;margin-bottom:14px;
  box-shadow:0 1px 4px rgba(0,0,0,0.06);}
.filter-panel h4{margin:0 0 10px;font-size:13px;color:__BLUE__;}
.filter-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:14px;}
.filter-item{display:flex;flex-direction:column;font-size:11px;color:__DARK__;}
.filter-item label{display:flex;justify-content:space-between;}
.filter-item input[type=range]{width:100%;accent-color:__BLUE__;}
.filter-counter{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:12px;margin-top:14px;}
.filter-counter .box{background:__BG__;border-radius:6px;padding:10px 12px;}
.filter-counter .big{font-size:22px;font-weight:700;color:__BLUE__;}
.filter-counter .label{font-size:10px;color:__GREY__;text-transform:uppercase;letter-spacing:0.04em;}
.filter-counter .means{font-size:11px;color:__GREY__;margin-top:6px;}
.filter-reset{margin-top:10px;font-size:11px;color:__BLUE__;cursor:pointer;background:none;border:1px solid __BLUE__;
  border-radius:14px;padding:3px 10px;}
.filter-reset:hover{background:__BLUE__;color:#fff;}
.filter-hint{font-size:10px;color:__GREY__;margin-top:2px;}
.filter-hint b{color:__ORANGE__;}
</style>
<div class="filter-panel">
  <h4>Filter samples — move sliders to live-update the QC pill in each row</h4>
  <div class="filter-grid">
    <div class="filter-item">
      <label>Sylph ANI ≥ <b id="lbl_ani">97 %</b></label>
      <input id="sl_ani" type="range" min="80" max="100" step="0.1" value="97">
      <div class="filter-hint">dataset min: <b id="hint_ani">—</b></div>
    </div>
    <div class="filter-item">
      <label>Sylph %cov ≥ <b id="lbl_cov">85 %</b></label>
      <input id="sl_cov" type="range" min="0" max="100" step="0.5" value="85">
      <div class="filter-hint">dataset min: <b id="hint_cov">—</b></div>
    </div>
    <div class="filter-item">
      <label>%alignment ≥ <b id="lbl_align">90 %</b></label>
      <input id="sl_align" type="range" min="0" max="100" step="0.5" value="90">
      <div class="filter-hint">dataset min: <b id="hint_align">—</b></div>
    </div>
    <div class="filter-item">
      <label>depth estimated ≥ <b id="lbl_depth">30×</b></label>
      <input id="sl_depth" type="range" min="0" max="__MAX_DEPTH__" step="1" value="30">
      <div class="filter-hint">dataset min: <b id="hint_depth">—</b></div>
    </div>
  </div>
  <div class="filter-counter">
    <div class="box"><div class="label">Total samples</div><div class="big" id="ct_vis">— / —</div></div>
    <div class="box"><div class="label">PASS my filter</div><div class="big" style="color:__GREEN__;" id="ct_pass">0</div></div>
    <div class="box"><div class="label">FAIL my filter</div><div class="big" style="color:__RED__;" id="ct_fail">0</div></div>
    <div class="box"><div class="label">Means (PASS)</div><div class="means" id="ct_means">—</div></div>
  </div>
  <button class="filter-reset" id="filter_reset">Reset to presets</button>
</div>
<script>
/* The filter panel sits ABOVE the table in source order, so the script body
   has to wait for DOMContentLoaded before it can find #qc_table.  Without
   this guard, document.getElementById('qc_table') returns null and update()
   throws, leaving all counters stuck at their placeholder ("— / —"). */
document.addEventListener('DOMContentLoaded', function(){
  var sl = {
    ani:   document.getElementById('sl_ani'),
    cov:   document.getElementById('sl_cov'),
    align: document.getElementById('sl_align'),
    depth: document.getElementById('sl_depth'),
  };
  var lbl = {
    ani:   document.getElementById('lbl_ani'),
    cov:   document.getElementById('lbl_cov'),
    align: document.getElementById('lbl_align'),
    depth: document.getElementById('lbl_depth'),
  };
  var ctVis  = document.getElementById('ct_vis');
  var ctPass = document.getElementById('ct_pass');
  var ctFail = document.getElementById('ct_fail');
  var ctMeans= document.getElementById('ct_means');
  var table  = document.getElementById('qc_table');
  var allRows = Array.from(table.tBodies[0].rows);
  var TOTAL = allRows.length;
  /* Preset values — what the "Reset" button restores them to. */
  var PRESET = { ani: 97, cov: 85, align: 90, depth: 30 };
  /* Compute and display dataset minimums per column so the user knows from
     where the slider starts having effect. */
  function minVal(key){
    var m = Infinity;
    allRows.forEach(function(r){ var v=+r.dataset[key]; if (!isNaN(v) && v<m) m=v; });
    return m === Infinity ? 0 : m;
  }
  document.getElementById('hint_ani').textContent   = minVal('ani').toFixed(2)   + ' %';
  document.getElementById('hint_cov').textContent   = minVal('cov').toFixed(1)   + ' %';
  document.getElementById('hint_align').textContent = minVal('align').toFixed(1) + ' %';
  document.getElementById('hint_depth').textContent = minVal('depth').toFixed(0) + '×';

  function update(){
    var ti = +sl.ani.value, tc = +sl.cov.value, ta = +sl.align.value, td = +sl.depth.value;
    lbl.ani.textContent   = ti.toFixed(2).replace(/\.?0+$/,'') + ' %';
    lbl.cov.textContent   = tc.toFixed(1).replace(/\.0$/,'') + ' %';
    lbl.align.textContent = ta.toFixed(1).replace(/\.0$/,'') + ' %';
    lbl.depth.textContent = td + '×';
    var pass=0, fail=0;
    var sI=0,sC=0,sA=0,sD=0;
    allRows.forEach(function(r){
      var ani = +r.dataset.ani,
          cov = +r.dataset.cov,
          a   = +r.dataset.align,
          d   = +r.dataset.depth;
      var ok = ani>=ti && cov>=tc && a>=ta && d>=td;
      /* Live-update the QC pill in this row instead of hiding the row.
         All 60 rows always remain visible; user sees how their slider
         choice flips each sample between PASS and FAIL in real time. */
      var pill = r.querySelector('.qc_pill');
      if (pill) {
        pill.classList.toggle('pass', ok);
        pill.classList.toggle('fail', !ok);
        pill.textContent = ok ? 'PASS' : 'FAIL';
      }
      if (ok) { pass++; sI+=ani; sC+=cov; sA+=a; sD+=d; }
      else    { fail++; }
    });
    ctVis.textContent  = TOTAL + ' / ' + TOTAL;
    ctPass.textContent = pass + ' (' + (100*pass/TOTAL).toFixed(1) + '%)';
    ctFail.textContent = fail;
    ctMeans.innerHTML = pass ?
      ('ANI <b>'+(sI/pass).toFixed(2)+'%</b> · cov <b>'+(sC/pass).toFixed(1)+'%</b><br>' +
       'align <b>'+(sA/pass).toFixed(1)+'%</b> · depth <b>'+(sD/pass).toFixed(0)+'×</b>')
      : '—';
  }
  Object.values(sl).forEach(function(s){ s.addEventListener('input', update); });
  document.getElementById('filter_reset').addEventListener('click', function(){
    sl.ani.value   = PRESET.ani;
    sl.cov.value   = PRESET.cov;
    sl.align.value = PRESET.align;
    sl.depth.value = PRESET.depth;
    update();
  });
  update();
});
</script>
"""


def _filter_panel_html(rows: list[SampleRow]) -> str:
    """Renders the slider + counter panel that controls qc_table visibility."""
    max_depth = max((r.est_depth for r in rows), default=0.0)
    return (FILTER_JS_TEMPLATE
            .replace("__BLUE__",  BLUE)
            .replace("__ORANGE__", ORANGE)
            .replace("__LBLUE__", LBLUE)
            .replace("__DARK__",  DARK)
            .replace("__BG__",    BG)
            .replace("__GREY__",  GREY)
            .replace("__GREEN__", GREEN)
            .replace("__RED__",   RED)
            .replace("__MAX_DEPTH__", str(int(max(50, np.ceil(max_depth))))))


def render_secondary_qc_html(rows: list[SampleRow], out: Path) -> None:
    """Standalone per-sample QC table with bar-fill cells, sortable columns,
    live filter sliders and a means/PASS/FAIL counter."""
    # Column min/max needed for the bar fills.
    def _rng(attr_or_callable):
        vals = []
        for r in rows:
            v = attr_or_callable(r) if callable(attr_or_callable) else getattr(r, attr_or_callable)
            if v: vals.append(float(v))
        return (min(vals) if vals else 0.0, max(vals) if vals else 0.0)

    rng_raw_reads = _rng("raw_illu_reads")
    rng_raw_gb    = _rng(lambda r: r.raw_gb)
    rng_clean_reads = _rng("clean_illu_reads")
    rng_clean_gb    = _rng(lambda r: r.clean_gb)
    rng_retention   = _rng(lambda r: r.retention_pct)
    rng_ani     = _rng("sylph_ani")
    rng_cov     = _rng("sylph_cov")
    rng_align   = _rng("align_rate")
    rng_depth   = _rng("est_depth")
    rng_contam  = _rng("contam")

    head = ("<table id='qc_table' class='dat sortable'><thead><tr>"
            "<th>Sample</th><th>Platform</th><th>Batch</th>"
            "<th>Raw reads</th><th>Raw GB</th>"
            "<th>Clean reads</th><th>Clean GB</th><th>%retention</th>"
            "<th>match_sylph</th>"
            "<th>Sylph ANI</th>"
            "<th>%cov_sylph</th><th>%alignment</th>"
            "<th>depth_estimated</th><th>%contam</th>"
            "<th>QC</th></tr></thead><tbody>")
    rows_html = []
    for r in rows:
        # Pill carries `qc_pill` so the filter JS can update it live as the
        # sliders move (no hiding of rows).  The td gets `title=` with the
        # pipeline's own QC verdict so the user can hover for the original.
        pill = (f"<span class='pill qc_pill {'pass' if r.qc_flag == 'PASS' else 'fail'}'>"
                f"{r.qc_flag}</span>")
        rows_html.append(
            "<tr"
            f" data-ani='{r.sylph_ani:.2f}'"
            f" data-cov='{r.sylph_cov:.1f}'"
            f" data-align='{r.align_rate:.1f}'"
            f" data-depth='{r.est_depth:.1f}'"
            f" data-retention='{r.retention_pct:.1f}'"
            f" data-clean-gb='{r.clean_gb:.3f}'"
            f" data-qc='{r.qc_flag}'>"
            f"<td>{r.sample_id}</td><td>{r.platform}</td><td>{r.batch or '—'}</td>"
            + _bar_cell(f"{r.raw_illu_reads:,}", r.raw_illu_reads, *rng_raw_reads, COLOUR_FASTQ)
            + _bar_cell(f"{r.raw_gb:.2f}",       r.raw_gb,         *rng_raw_gb,    COLOUR_FASTQ)
            + _bar_cell(f"{r.clean_illu_reads:,}", r.clean_illu_reads, *rng_clean_reads, COLOUR_FASTQ)
            + _bar_cell(f"{r.clean_gb:.2f}",     r.clean_gb,       *rng_clean_gb,  COLOUR_FASTQ)
            + _bar_cell(f"{r.retention_pct:.1f}", r.retention_pct, *rng_retention, COLOUR_FASTQ)
            + f"<td>{r.species or '—'}</td>"
            + _bar_cell(f"{r.sylph_ani:.2f}",  r.sylph_ani,  *rng_ani,   COLOUR_SYLPH)
            + _bar_cell(f"{r.sylph_cov:.1f}",  r.sylph_cov,  *rng_cov,   COLOUR_SYLPH)
            + _bar_cell(f"{r.align_rate:.1f}", r.align_rate, *rng_align, COLOUR_SYLPH)
            + _bar_cell(f"{r.est_depth:.1f}",  r.est_depth,  *rng_depth, COLOUR_SYLPH)
            + _bar_cell(f"{r.contam:.2f}",     r.contam,     *rng_contam, COLOUR_SYLPH)
            + f"<td class='qc_cell' title='Pipeline QC: {r.qc_flag}'>{pill}</td></tr>")
    body = head + "".join(rows_html) + "</tbody></table>"

    html = (f"<!DOCTYPE html><html><head><meta charset='utf-8'>"
            f"<title>Fastq QC — per sample</title><style>{CSS}</style></head>"
            f"<body><div class='wrap'><section class='section'>"
            f"<h3>Fastq QC — per-sample table</h3>"
            f"<p class='muted'>Bar fill inside each numeric cell is proportional "
            f"to the value within its column. Blue = fastq metrics, orange = identification &amp; coverage. "
            f"Click any column header to sort.</p>"
            f"{_filter_panel_html(rows)}"
            f"{body}"
            f"{SORT_JS}"
            f"</section></div></body></html>")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")


# ─── Phase 4 — BUSCO + QUAST ───────────────────────────────────────────────

# All metrics here are restricted to contigs ≥ 1 kb (the standard size cutoff
# for short-read assemblies — sub-1 kb fragments are usually mis-assemblies
# or low-coverage debris and dilute the headline numbers).
QUAST_KEY_COLS = [
    ("Total length (>= 1000 bp)",  "Total length"),
    ("# contigs (>= 1000 bp)",     "# Contigs"),
    ("Largest contig",             "Largest contig"),
    ("GC (%)",                     "GC%"),
    ("N50",                        "N50"),
    ("L50",                        "L50"),
    ("Genome fraction (%)",        "Genome fraction (%)"),
    ("Misassembled contigs length","Misassembled length"),
    ("# mismatches per 100 kbp",   "Mismatches/100Kb"),
    ("Duplication ratio",          "Duplication ratio"),
]

BUSCO_RE = re.compile(
    r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%")


def collect_quast(results_dir: Path) -> list[dict]:
    out = []
    for tsv in sorted((results_dir / "05_post_asm_qc").glob("*/quast/*.quast/transposed_report.tsv")):
        recs = _read_tsv(tsv)
        if not recs: continue
        rec = recs[0]
        sid = rec.get("Assembly", tsv.parent.parent.parent.name)
        d = {"Sample": sid}
        for k, label in QUAST_KEY_COLS:
            d[label] = rec.get(k, "—")
        out.append(d)
    return out


def collect_busco(results_dir: Path) -> list[dict]:
    out = []
    for txt in sorted((results_dir / "05_post_asm_qc").glob("*/busco/short_summary.specific*.txt")):
        sid = txt.parent.parent.name
        content = txt.read_text(encoding="utf-8", errors="ignore")
        m = BUSCO_RE.search(content)
        if not m: continue
        out.append({
            "Sample": sid,
            "C":      float(m.group(1)),
            "S":      float(m.group(2)),
            "D":      float(m.group(3)),
            "F":      float(m.group(4)),
            "M":      float(m.group(5)),
        })
    return out


def _quast_numeric(value: str) -> Optional[float]:
    """Parse QUAST cell text into a float; returns None when not numeric."""
    if value is None: return None
    v = str(value).strip()
    if v in ("", "—", "-", "NA", "nan"): return None
    # Many QUAST values include compound forms like '1 + 4 part'; keep the first
    # numeric token so the cell can still feed the bar-fill rendering.
    try:
        return float(v.split()[0])
    except (TypeError, ValueError):
        return None


def render_assembly_qc(quast_rows: list[dict], busco_rows: list[dict],
                       secondary_html_path: Optional[Path]) -> str:
    """Main HTML assembly section.

    Shows aggregate (mean / median) assembly stats + a Plotly BUSCO chart
    with per-sample hover.  The full QUAST table is offloaded to a secondary
    HTML linked from the top of the section.
    """
    if not quast_rows and not busco_rows:
        return ""

    # ── Aggregate stats card.
    def _agg(col):
        vs = [_quast_numeric(r.get(col)) for r in quast_rows]
        vs = [v for v in vs if v is not None]
        return (np.mean(vs) if vs else 0.0, np.median(vs) if vs else 0.0,
                np.min(vs) if vs else 0.0,  np.max(vs) if vs else 0.0)

    fmt_size  = lambda v: f"{v/1e6:.2f} Mb" if v else "—"
    fmt_int   = lambda v: f"{int(v):,}"     if v else "—"
    fmt_pct   = lambda v: f"{v:.1f} %"      if v else "—"

    def stat_card(label, mean, median, mn, mx, fmt):
        return (f'<div class="card"><div class="label">{label}</div>'
                f'<div class="big" style="font-size:18px;">{fmt(mean)}</div>'
                f'<div class="muted" style="font-size:11px;">'
                f'median {fmt(median)} &middot; range {fmt(mn)}–{fmt(mx)}</div></div>')

    cards = []
    if quast_rows:
        mean, med, mn, mx = _agg("Total length")
        cards.append(stat_card("Total length",         mean, med, mn, mx, fmt_size))
        mean, med, mn, mx = _agg("# Contigs")
        cards.append(stat_card("# Contigs",            mean, med, mn, mx, fmt_int))
        # Mean contig length per sample, then aggregated across the cohort.
        mean_lens = []
        for r in quast_rows:
            tl = _quast_numeric(r.get("Total length"))
            nc = _quast_numeric(r.get("# Contigs"))
            if tl and nc and nc > 0:
                mean_lens.append(tl / nc)
        if mean_lens:
            cards.append(stat_card("Mean contig length",
                                   float(np.mean(mean_lens)),
                                   float(np.median(mean_lens)),
                                   float(np.min(mean_lens)),
                                   float(np.max(mean_lens)), fmt_size))
        mean, med, mn, mx = _agg("Largest contig")
        cards.append(stat_card("Largest contig",       mean, med, mn, mx, fmt_size))
        mean, med, mn, mx = _agg("Genome fraction (%)")
        cards.append(stat_card("Genome fraction",      mean, med, mn, mx, fmt_pct))
        mean, med, mn, mx = _agg("GC%")
        cards.append(stat_card("GC content",           mean, med, mn, mx, fmt_pct))
        mean, med, mn, mx = _agg("Mismatches/100Kb")
        cards.append(stat_card("Mismatches / 100 Kb",  mean, med, mn, mx, lambda v: f"{v:.1f}"))
    if busco_rows:
        Cs = [r["C"] for r in busco_rows]
        cards.append(stat_card("BUSCO complete (C)",
                               float(np.mean(Cs)), float(np.median(Cs)),
                               float(np.min(Cs)),  float(np.max(Cs)), fmt_pct))
    stats_html = (f'<p class="muted" style="margin-bottom:6px;">All assembly '
                  f'metrics computed on contigs ≥ 1 kb. Cards show '
                  f'<b>mean</b> with median &middot; range underneath.</p>'
                  f'<div class="cards">{"".join(cards)}</div>')

    # ── Plotly BUSCO chart with per-sample hover.
    busco_chart = ""
    if busco_rows:
        b = sorted(busco_rows, key=lambda x: x["Sample"])
        sids = [r["Sample"] for r in b]
        traces = []
        for key, colour, label in [("S", GREEN,  "Complete-single"),
                                   ("D", BLUE,   "Complete-duplicated"),
                                   ("F", ORANGE, "Fragmented"),
                                   ("M", RED,    "Missing")]:
            ys = [r[key] for r in b]
            traces.append({
                "type": "bar", "x": sids, "y": ys, "name": label,
                "marker": {"color": colour, "line": {"color": "white", "width": 0.3}},
                "hovertemplate": f"<b>%{{x}}</b><br>{label}: %{{y:.2f}}%<extra></extra>",
            })
        layout = {
            "barmode": "stack",
            "title": {"text": "BUSCO completeness (saccharomycetes_odb12)",
                      "font": {"size": 14, "color": BLUE}},
            "xaxis": {"title": "Sample", "tickangle": -90,
                      "automargin": True, "tickfont": {"size": 9}},
            "yaxis": {"title": "% BUSCOs", "range": [0, 101]},
            "legend": {"orientation": "h", "x": 0, "y": -0.3},
            "margin": {"l": 60, "r": 30, "t": 60, "b": 140},
            "plot_bgcolor": "white", "paper_bgcolor": "white",
        }
        busco_chart = _plotly_div(traces, layout,
                                  height_px=max(420, 0.18 * 6 * len(sids) + 300))

    link_html = ""
    if secondary_html_path:
        link_html = (f'<p><a href="{secondary_html_path.name}" target="_blank" '
                     f'style="color:{BLUE};font-weight:600;">'
                     f'Open full per-sample assembly QC table →</a></p>')

    return f"""
    <section class="section">
      <h3>Assembly QC — BUSCO &amp; QUAST</h3>
      {link_html}
      <h4>Aggregate metrics</h4>
      {stats_html}
      <h4>BUSCO completeness</h4>
      {busco_chart}
    </section>
    """


def render_secondary_assembly_html(quast_rows: list[dict],
                                   busco_rows: list[dict],
                                   out: Path) -> None:
    """Standalone HTML with the full QUAST + BUSCO per-sample table.
    Numeric cells are shaded by bar fill so rows can be compared at a glance.
    """
    if not quast_rows and not busco_rows:
        return

    # Column ranges across the cohort for the bar-fill normalisation.
    quast_rng = {label: (
        (lambda vs: (min(vs), max(vs)))(
            [v for v in (_quast_numeric(r.get(label)) for r in quast_rows) if v is not None]
            or [0.0, 0.0]))
        for _, label in QUAST_KEY_COLS}

    busco_by_sample = {r["Sample"]: r for r in busco_rows}
    busco_keys = [("C", "BUSCO C %"), ("S", "BUSCO S %"),
                  ("D", "BUSCO D %"), ("F", "BUSCO F %"), ("M", "BUSCO M %")]
    busco_rng = {k: (0.0, 100.0) for k, _ in busco_keys}

    head = "<table class='dat sortable'><thead><tr><th>Sample</th>"
    for _, lbl in QUAST_KEY_COLS:
        head += f"<th>{lbl}</th>"
    for _, lbl in busco_keys:
        head += f"<th>{lbl}</th>"
    head += "</tr></thead><tbody>"

    body_rows = []
    for r in quast_rows:
        cells = [f"<td>{r['Sample']}</td>"]
        for raw_col, lbl in QUAST_KEY_COLS:
            txt = r.get(lbl, "—")
            num = _quast_numeric(txt)
            if num is None:
                cells.append(f"<td>{txt}</td>")
            else:
                vmin, vmax = quast_rng[lbl]
                cells.append(_bar_cell(str(txt), num, vmin, vmax, COLOUR_ASSEMBLY))
        bu = busco_by_sample.get(r["Sample"], {})
        for k, _ in busco_keys:
            v = bu.get(k)
            if v is None:
                cells.append("<td>—</td>")
            else:
                cells.append(_bar_cell(f"{v:.1f}", float(v), *busco_rng[k], COLOUR_ASSEMBLY))
        body_rows.append("<tr>" + "".join(cells) + "</tr>")
    body = head + "".join(body_rows) + "</tbody></table>"

    html = (f"<!DOCTYPE html><html><head><meta charset='utf-8'>"
            f"<title>Assembly QC — per sample</title><style>{CSS}</style></head>"
            f"<body><div class='wrap'><section class='section'>"
            f"<h3>Assembly QC — per-sample table</h3>"
            f"<p class='muted'>Bar fill in each numeric cell is proportional to "
            f"its value within the column (green = assembly / QC metrics). "
            f"Click any column header to sort.</p>"
            f"<div style='overflow-x:auto;'>{body}</div>"
            f"{SORT_JS}"
            f"</section></div></body></html>")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")


# ─── Phase 5 — AMR + CNV heatmaps (embed existing PNGs) ────────────────────

# ─── CNV log2 heatmap (in-process, mirrors cnvkit_heatmap.py template) ─────
#
# The pipeline's CNV_VISUALIZATION module sometimes fails to publish PNGs.
# The report keeps a fallback that reads the aggregated log2 matrix plus the
# per-sample cnvkit_cn from cnv_events.tsv and renders the divergent heatmap.

def _ploidy_from_results(results_dir: Path) -> int:
    """Most common ploidy across species_call.tsv; defaults to 1 (C. auris)."""
    counts: dict[int, int] = {}
    for sp in (results_dir / "03_identification").glob("*/*.species_call.tsv"):
        recs = _read_tsv(sp)
        if not recs: continue
        v = recs[0].get("ploidy")
        try: p = int(float(v)) if v else 1
        except (TypeError, ValueError): p = 1
        counts[p] = counts.get(p, 0) + 1
    if not counts: return 1
    return max(counts.items(), key=lambda kv: kv[1])[0]


def render_cnv_heatmap_inline(results_dir: Path) -> str:
    """Divergent log2 heatmap (samples × genes) with cnvkit_cn annotations.
    Reuses the cohort log2 matrix + per-sample cnv_events.tsv files.
    """
    try:
        import pandas as pd
        from matplotlib.colors import TwoSlopeNorm
    except ImportError:
        return ""

    log2_path = results_dir / "08_cnv" / "aggregated" / "cnv_log2_matrix.tsv"
    if not log2_path.exists():
        return ""
    log2 = pd.read_csv(log2_path, sep="\t").set_index("sample_id")
    if log2.empty:
        return ""

    # Per-sample cn + event matrices (sparse — only filled where we have an event).
    cn  = pd.DataFrame(np.nan, index=log2.index, columns=log2.columns, dtype=float)
    evm = pd.DataFrame("",   index=log2.index, columns=log2.columns, dtype=object)
    for sid in log2.index:
        ev_path = results_dir / "08_cnv" / sid / f"{sid}.cnv_events.tsv"
        if not ev_path.exists(): continue
        try:
            ev = pd.read_csv(ev_path, sep="\t")
        except Exception:
            continue
        if "gene" not in ev.columns: continue
        for _, row in ev.iterrows():
            g = str(row.get("gene", "")).strip()
            if g in evm.columns:
                evm.loc[sid, g] = str(row.get("event", "") or "").strip().lower()
            v = row.get("cnvkit_cn", row.get("cn"))
            if g in cn.columns and v is not None and str(v) not in ("", "nan", "NA"):
                try: cn.loc[sid, g] = float(v)
                except (TypeError, ValueError): pass

    # Drop gene columns where no sample has a non-`normal` event.  When the
    # cohort has zero real events across every gene, skip the panel entirely
    # — a heatmap of "all-normal" is not informative and just adds noise.
    EMPTY = {"", "normal", "nan", "no_cov", "none"}
    keep_genes = [g for g in log2.columns
                  if not evm[g].astype(str).isin(EMPTY).all()]
    if not keep_genes:
        return ""
    if len(keep_genes) < len(log2.columns):
        log2 = log2.loc[:, keep_genes]
        cn   = cn.loc[:, keep_genes]
        evm  = evm.loc[:, keep_genes]

    ploidy = _ploidy_from_results(results_dir)
    M = log2.values.astype(float)
    vmax = max(1.0, float(np.nanmax(np.abs(M)))) if M.size else 1.0
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    n_samples, n_genes = M.shape
    fig_w = max(6, 0.55 * n_genes + 2.5)
    fig_h = max(4, 0.12 * n_samples + 1.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(M, cmap="RdBu_r", norm=norm, aspect="auto",
                   interpolation="nearest")

    # Annotate cells where cn departs from ploidy (cleaner than annotating all).
    for i in range(n_samples):
        for j in range(n_genes):
            v_log2 = M[i, j]
            v_cn   = cn.iloc[i, j]
            if pd.isna(v_cn): continue
            cn_i  = int(round(float(v_cn)))   # round, don't truncate
            delta = cn_i - ploidy
            # Skip near-zero quiet cells.
            if delta == 0 and abs(v_log2) < 0.30: continue
            label = (f"{cn_i}" if delta == 0
                     else f"{cn_i}({'+' if delta > 0 else '−'}{abs(delta)})")
            txt_col = "white" if abs(v_log2) > vmax * 0.55 else DARK
            ax.text(j, i, label, ha="center", va="center",
                    fontsize=5, color=txt_col)

    ax.set_xticks(range(n_genes))
    ax.set_xticklabels(list(log2.columns), rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(n_samples))
    ax.set_yticklabels(list(log2.index), fontsize=5)
    ax.set_xlabel(""); ax.set_ylabel("Sample")
    cbar = plt.colorbar(im, ax=ax, label="CNVkit log2 ratio",
                        shrink=0.5, pad=0.02)
    cbar.ax.tick_params(labelsize=9)
    ax.set_title(f"CNV — CNVkit log2 ratio (n={n_samples} samples, ploidy {ploidy})",
                 fontsize=12, color=BLUE)
    fig.tight_layout()

    legend = (
        '<div style="display:flex;gap:18px;margin:6px 0 4px;font-size:11px;color:'
        f'{DARK};">'
        f'<span><span style="display:inline-block;width:14px;height:14px;background:#b2182b;border-radius:3px;vertical-align:middle;margin-right:5px;"></span>gain (log2 > 0)</span>'
        f'<span><span style="display:inline-block;width:14px;height:14px;background:#f7f7f7;border:1px solid #ccc;border-radius:3px;vertical-align:middle;margin-right:5px;"></span>neutral</span>'
        f'<span><span style="display:inline-block;width:14px;height:14px;background:#2166ac;border-radius:3px;vertical-align:middle;margin-right:5px;"></span>loss (log2 &lt; 0)</span>'
        '</div>')
    return ("<h4>CNV — CNVkit log2 heatmap</h4>" + legend
            + f'<img class="fig" src="{_fig_to_b64(fig, dpi=140)}" alt="cnv_log2">')


def render_amr_cnv_combined_inline(results_dir: Path,
                                   rows: list[SampleRow],
                                   batch_colors: dict[str, str]) -> str:
    """Single dual-panel heatmap: AMR (left) + CNV log2 (right) sharing the
    sample axis.  Built inline so we can guarantee:
      · R / I / S calls overlaid as letters on each cell (so 'I' is visible),
      · sample rows grouped by batch with a batch-coloured strip + Y labels,
      · a clear top legend covering both panels.
    """
    try:
        import pandas as pd
        from matplotlib.colors import ListedColormap, TwoSlopeNorm
        from matplotlib.gridspec import GridSpec
        import matplotlib.patches as mpatches
    except ImportError:
        return ""

    amr_path = results_dir / "cohort" / "call_matrix.tsv"
    cnv_path = results_dir / "08_cnv" / "aggregated" / "cnv_log2_matrix.tsv"
    if not amr_path.exists() or not cnv_path.exists():
        return ""
    amr = pd.read_csv(amr_path, sep="\t").set_index("sample_id")
    cnv = pd.read_csv(cnv_path, sep="\t").set_index("sample_id")
    common = sorted(set(amr.index) & set(cnv.index))
    if not common:
        return ""

    # Order samples by batch (alphabetical) then sample id.
    sample_batch = {r.sample_id: r.batch for r in rows}
    ordered = sorted(common, key=lambda s: ((sample_batch.get(s) or "~zzz"), s))
    amr = amr.loc[ordered]
    cnv = cnv.loc[ordered]

    n = len(ordered)
    n_drugs = amr.shape[1]
    n_genes = cnv.shape[1]
    call_to_num = {"R": 2, "I": 1, "S": 0}
    # Vectorised numeric encoding (compatible with pandas <2.1, no DataFrame.map).
    amr_num = pd.DataFrame(
        np.vectorize(lambda v: call_to_num.get(str(v).strip(), 0))(amr.values),
        index=amr.index, columns=amr.columns)

    ploidy = _ploidy_from_results(results_dir)

    # CNV cn matrix (for cell labels) + event matrix (for color masking).
    # The pipeline's event column ('normal', 'dup_low', 'del', ...) is the
    # authoritative call: we suppress colour on cells flagged 'normal' so
    # that mosdepth log2 noise (typically ±0.5 around 0) does NOT appear as
    # a fake loss/gain on the heatmap.
    cn  = pd.DataFrame(np.nan, index=ordered, columns=cnv.columns, dtype=float)
    evm = pd.DataFrame("",   index=ordered, columns=cnv.columns, dtype=object)
    for sid in ordered:
        ev = results_dir / "08_cnv" / sid / f"{sid}.cnv_events.tsv"
        if not ev.exists(): continue
        try:
            df = pd.read_csv(ev, sep="\t")
        except Exception:
            continue
        for _, row in df.iterrows():
            g = str(row.get("gene", "")).strip()
            event = str(row.get("event", "") or "").strip().lower()
            if g in evm.columns:
                evm.loc[sid, g] = event
            # Prefer cnvkit_cn / cn directly; otherwise estimate from norm × ploidy
            # (the cnv_events.tsv emitted by the pipeline only ships `norm` and
            # `log2_ratio`, so a direct cn column is usually absent).
            v = row.get("cnvkit_cn", row.get("cn"))
            if v is None or str(v).strip() in ("", "nan", "NA"):
                norm_val = row.get("norm")
                if norm_val is not None and str(norm_val).strip() not in ("", "nan", "NA"):
                    try: v = float(norm_val) * ploidy
                    except (TypeError, ValueError): v = None
            if g in cn.columns and v is not None and str(v) not in ("", "nan", "NA"):
                try: cn.loc[sid, g] = float(v)
                except (TypeError, ValueError): pass

    # Build the colour matrix: zero out cells whose event is 'normal' (or
    # blank) so they paint pure white, regardless of mosdepth log2 noise.
    M = cnv.values.astype(float).copy()
    for i, sid in enumerate(ordered):
        for j, g in enumerate(cnv.columns):
            ev_str = evm.loc[sid, g]
            if ev_str in ("", "normal", "nan"):
                M[i, j] = 0.0

    # Drop gene columns where *no* sample has a real (≠ normal) event so the
    # CNV panel stays compact.  When a new gene gains an event in any later
    # cohort it shows up automatically; until then the column is hidden.
    # If the cohort has zero real events across every gene the combined
    # heatmap is skipped (we'd just render an all-normal CNV panel next to
    # AMR, which is misleading).
    EMPTY = {"", "normal", "nan", "no_cov", "none"}
    keep_genes = [g for g in cnv.columns
                  if not evm[g].astype(str).str.strip().str.lower().isin(EMPTY).all()]
    if not keep_genes:
        return ""
    if len(keep_genes) < len(cnv.columns):
        keep_idx = [list(cnv.columns).index(g) for g in keep_genes]
        cnv = cnv.loc[:, keep_genes]
        cn  = cn.loc[:, keep_genes]
        evm = evm.loc[:, keep_genes]
        M   = M[:, keep_idx]
        n_genes = len(keep_genes)

    vmax = max(1.0, float(np.nanmax(np.abs(M)))) if M.size else 1.0
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    # Wider strip column so that full sample IDs (~9-12 chars) fit.
    strip_units = 2.5
    fig_w = max(13, 0.45 * (n_drugs + n_genes) + strip_units + 4.5)
    fig_h = max(7,  0.16 * n + 2.5)
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = GridSpec(1, 4, figure=fig,
                  width_ratios=[strip_units, n_drugs, n_genes, 0.55],
                  wspace=0.02, left=0.04, right=0.94, top=0.92, bottom=0.18)
    ax_b = fig.add_subplot(gs[0, 0])
    ax_a = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])
    ax_cbar = fig.add_subplot(gs[0, 3])

    # Batch strip + sample names — text larger and centred over the wider strip.
    label_font = max(7, min(11, 900 // n))
    for i, s in enumerate(ordered):
        b = sample_batch.get(s, "")
        col = _hex_for_batch(b or "n/a", batch_colors) if b else "#dddddd"
        ax_b.add_patch(mpatches.Rectangle(
            (-0.5, i - 0.5), 1, 1, facecolor=col, edgecolor="none"))
        lum = _hex_to_lum(col)
        ax_b.text(0, i, s, fontsize=label_font, ha="center", va="center",
                  color="white" if lum < 0.55 else "#1a1a1a")
    ax_b.set_xlim(-0.5, 0.5); ax_b.set_ylim(n - 0.5, -0.5)
    ax_b.set_xticks([]); ax_b.set_yticks([])
    for sp in ax_b.spines.values(): sp.set_visible(False)

    # AMR panel — discrete colormap S / I / R + letter overlay on each cell.
    AMR_CMAP = ListedColormap([GREEN, ORANGE, RED])
    ax_a.imshow(amr_num.values, cmap=AMR_CMAP, vmin=0, vmax=2,
                aspect="auto", interpolation="nearest")
    cell_font = max(5, min(9, 600 // n))
    for i in range(n):
        for j in range(n_drugs):
            v = str(amr.iloc[i, j]).strip()
            if v in ("R", "I"):
                ax_a.text(j, i, v, ha="center", va="center",
                          fontsize=cell_font, color="white", weight="bold")
    ax_a.set_xticks(range(n_drugs))
    ax_a.set_xticklabels(amr.columns, rotation=45, ha="right", fontsize=12,
                         fontweight="medium")
    ax_a.set_yticks([])
    ax_a.set_title("AMR (R / I / S)", fontsize=11, color=BLUE, pad=4)

    # CNV panel — divergent + cn delta annotations.
    im = ax_c.imshow(M, cmap="RdBu_r", norm=norm, aspect="auto",
                     interpolation="nearest")
    for i in range(n):
        for j in range(n_genes):
            v_cn = cn.iloc[i, j]
            if pd.isna(v_cn): continue
            # round() — int() truncates toward zero, so norm=0.97 (a normal
            # haploid cell) would become cn=0 and falsely label "0cn (−1)".
            cn_i = int(round(float(v_cn))); delta = cn_i - ploidy
            if delta == 0:
                # Normal cell — leave blank.  The heatmap colour (white /
                # near-white) already communicates "no change vs ploidy".
                continue
            sign = "+" if delta > 0 else "−"
            label = f"{cn_i}cn ({sign}{abs(delta)})"
            txt = "white" if abs(M[i, j]) > vmax * 0.55 else DARK
            ax_c.text(j, i, label, ha="center", va="center",
                      fontsize=cell_font, color=txt, weight="bold")
    ax_c.set_xticks(range(n_genes))
    ax_c.set_xticklabels(cnv.columns, rotation=45, ha="right", fontsize=12,
                         fontweight="medium")
    ax_c.set_yticks([])
    ax_c.set_title(f"CNV — CNVkit log2 (ploidy {ploidy})",
                   fontsize=11, color=BLUE, pad=4)

    plt.colorbar(im, cax=ax_cbar, label="CNV log2 ratio")
    fig.suptitle("AMR × CNV combined — samples × {drugs + genes}",
                 fontsize=13, weight="bold", color=BLUE, y=0.97)

    legend = (
        '<div style="display:flex;gap:18px;flex-wrap:wrap;margin:6px 0 4px;font-size:11px;color:'
        f'{DARK};">'
        f'<span><span style="display:inline-block;width:14px;height:14px;background:{RED};border-radius:3px;vertical-align:middle;margin-right:5px;"></span><b>R</b> · Resistant</span>'
        f'<span><span style="display:inline-block;width:14px;height:14px;background:{ORANGE};border-radius:3px;vertical-align:middle;margin-right:5px;"></span><b>I</b> · Intermediate</span>'
        f'<span><span style="display:inline-block;width:14px;height:14px;background:{GREEN};border-radius:3px;vertical-align:middle;margin-right:5px;"></span><b>S</b> · Susceptible</span>'
        '<span style="border-left:1px solid #ccc;padding-left:18px;margin-left:6px;">'
        '<span style="display:inline-block;width:14px;height:14px;background:#b2182b;border-radius:3px;vertical-align:middle;margin-right:5px;"></span>gain (log2 &gt; 0)</span>'
        '<span><span style="display:inline-block;width:14px;height:14px;background:#f7f7f7;border:1px solid #ccc;border-radius:3px;vertical-align:middle;margin-right:5px;"></span>neutral</span>'
        f'<span><span style="display:inline-block;width:14px;height:14px;background:#2166ac;border-radius:3px;vertical-align:middle;margin-right:5px;"></span>loss (log2 &lt; 0)</span>'
        '</div>')
    return (legend + f'<img class="fig" src="{_fig_to_b64(fig, dpi=140)}" alt="amr_cnv_combined">')


def render_amr_cnv(results_dir: Path, rows: list[SampleRow],
                   batch_colors: dict[str, str]) -> str:
    """Single combined AMR × CNV heatmap (per user request)."""
    combined = render_amr_cnv_combined_inline(results_dir, rows, batch_colors)
    if not combined:
        return ""
    return f"""
    <section class="section">
      <h3>AMR &amp; CNV</h3>
      <h4>AMR × CNV combined heatmap</h4>
      <p class="muted">Single dual-panel view. Left: AMR calls (S / I / R)
      with explicit letter overlay so intermediate is always visible. Right:
      CNVkit log2 ratio — <b>only cells with a real gain/loss carry a label</b>
      (format <code>Xcn (+/−N)</code>); blank/white cells = normal cn at the
      expected ploidy. Samples ordered by batch then alphabetical; batch
      strip on the left uses the same palette as the rest of the report.</p>
      {combined}
    </section>
    """


# ─── Phase 6 — SNP distance matrix + UPGMA / MSN network ───────────────────
#
# Inline implementations of the make_snp_heatmap.py and make_msn_upgma.py
# templates (issues/report/scripts_template/). Both are reproduced here so
# that the report module is self-contained and can colour the MSN by the same
# batch palette used elsewhere in the document (`batch_colors` registry).

def _load_snp_matrix(molten: Path) -> tuple[list[str], np.ndarray]:
    """Read a molten SNP distance file and return (samples_sorted, NxN ndarray)."""
    with open(molten) as fh:
        sep = "\t" if "\t" in fh.readline() else ","
    recs = []
    with open(molten) as fh:
        for line in fh:
            parts = line.rstrip().split(sep)
            if len(parts) < 3: continue
            try:
                recs.append((parts[0], parts[1], int(parts[2])))
            except ValueError:
                pass
    samples = sorted({a for a, _, _ in recs} | {b for _, b, _ in recs})
    idx = {s: i for i, s in enumerate(samples)}
    n = len(samples)
    M = np.zeros((n, n), dtype=int)
    for a, b, d in recs:
        M[idx[a], idx[b]] = d
    # Symmetrise (snpdists output is already symmetric but be defensive).
    M = np.maximum(M, M.T)
    return samples, M


def _hex_to_lum(h: str) -> float:
    h = h.lstrip("#")
    if len(h) != 6: return 0.5
    r, g, b = int(h[0:2], 16)/255, int(h[2:4], 16)/255, int(h[4:6], 16)/255
    return 0.299*r + 0.587*g + 0.114*b


def _render_snp_heatmap_plotly(samples: list[str], M: np.ndarray, cohort: str,
                               sample_batch_map: dict,
                               batch_colors: dict) -> str:
    """Plotly square SNP distance heatmap with UPGMA-driven sample order.

    - Hover on each cell shows the two samples + the integer SNP distance.
    - X-axis labels rotated 45°.
    - Y-axis labels only on the left (default Plotly).
    - Sample order: UPGMA leaves (average linkage on the SNP matrix).
    - Batch annotation strip drawn as a left margin column using shapes.
    """
    from scipy.cluster.hierarchy import linkage, dendrogram as scipy_dendr
    from scipy.spatial.distance import squareform

    n = len(samples)
    if n < 2: return ""
    Z = linkage(squareform(M.astype(float), checks=False), method="average")
    dn = scipy_dendr(Z, no_plot=True)
    order = dn["leaves"]
    samples_ord = [samples[i] for i in order]
    M_ord = M[np.ix_(order, order)]

    vmax = max(1, min(200, int(M.max())))
    # Yellow → garnet diverging scale (PowerNorm-like with manual stops).
    colorscale = [
        [0.00, "#fffde7"], [0.06, "#fff59d"], [0.14, "#fff176"],
        [0.25, "#ffd54f"], [0.40, "#ffa726"], [0.55, "#fb8c00"],
        [0.70, "#e53935"], [0.85, "#b71c1c"], [1.00, "#5e0e0e"],
    ]
    # Custom text labels for hover.
    customdata = [[samples_ord[i] + " ↔ " + samples_ord[j] for j in range(n)]
                  for i in range(n)]
    heatmap_trace = {
        "type": "heatmap",
        "z": M_ord.tolist(),
        "x": samples_ord, "y": samples_ord,
        "colorscale": colorscale, "zmin": 0, "zmax": vmax,
        "customdata": customdata,
        "hovertemplate": "%{customdata}<br><b>%{z} SNPs</b><extra></extra>",
        "colorbar": {"title": {"text": f"SNP distance<br>(cap {vmax})", "side": "right"},
                     "thickness": 12, "len": 0.55, "x": 1.02,
                     "tickfont": {"size": 10}},
        "showscale": True,
    }

    # Batch annotation overlaid as left strip via shapes outside the plot area.
    # We add coloured rectangles right next to the y-axis tick labels.
    shapes = []
    if any(sample_batch_map.get(s) for s in samples_ord):
        for i, s in enumerate(samples_ord):
            b = sample_batch_map.get(s, "")
            if not b: continue
            col = _hex_for_batch(b, batch_colors)
            shapes.append({
                "type": "rect", "xref": "paper", "yref": "y",
                "x0": -0.035, "x1": -0.015,
                "y0": i - 0.5, "y1": i + 0.5,
                "fillcolor": col, "line": {"width": 0},
            })

    cell_font = max(6, min(11, 1100 // n))
    layout = {
        "title": {"text": f"SNP distance heatmap — {cohort} (n={n}) · UPGMA-ordered",
                  "font": {"size": 13, "color": BLUE}},
        "xaxis": {"tickangle": -45, "tickfont": {"size": cell_font},
                  "side": "bottom", "automargin": True,
                  "showgrid": False, "ticks": ""},
        "yaxis": {"autorange": "reversed", "tickfont": {"size": cell_font},
                  "side": "left", "automargin": True, "showgrid": False,
                  "ticks": ""},
        "shapes": shapes,
        "margin": {"l": 130, "r": 100, "t": 60, "b": 130},
        "plot_bgcolor": "white", "paper_bgcolor": "white",
    }
    return _plotly_div([heatmap_trace], layout,
                       height_px=max(500, 14 * n + 200))


def _render_snp_heatmap_template(samples: list[str], M: np.ndarray, cohort: str,
                                 sample_batch_map: dict, batch_colors: dict) -> str:
    """Legacy matplotlib renderer — kept for fallback if Plotly is undesired.
    """
    from scipy.cluster.hierarchy import linkage, dendrogram as scipy_dendr
    from scipy.spatial.distance import squareform
    from matplotlib.colors import LinearSegmentedColormap, PowerNorm
    from matplotlib.gridspec import GridSpec
    import matplotlib.patches as mpatches

    n = len(samples)
    if n < 2: return ""

    Z = linkage(squareform(M.astype(float), checks=False), method="average")
    dn = scipy_dendr(Z, no_plot=True)
    order = dn["leaves"]
    samples_ord = [samples[i] for i in order]
    M_ord = M[np.ix_(order, order)]

    cell_size = max(0.20, min(0.45, 14.0 / n))
    fig_w = max(11, cell_size * n + 4.5)
    fig_h = max(8,  cell_size * n + 3.0)
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.0, 0.45, 8.0],
                  wspace=0.02, left=0.04, right=0.86, top=0.92, bottom=0.10)
    ax_dendro = fig.add_subplot(gs[0, 0])
    ax_annot  = fig.add_subplot(gs[0, 1])
    ax_heat   = fig.add_subplot(gs[0, 2])

    SNP_CMAP = LinearSegmentedColormap.from_list("snp_dist", [
        (0.00, "#fffde7"), (0.06, "#fff59d"), (0.14, "#fff176"),
        (0.25, "#ffd54f"), (0.40, "#ffa726"), (0.55, "#fb8c00"),
        (0.70, "#e53935"), (0.85, "#b71c1c"), (1.00, "#5e0e0e"),
    ])
    vmax = max(1, min(200, int(np.max(M))))
    norm = PowerNorm(gamma=0.55, vmin=0, vmax=vmax)

    # Left dendrogram (rotated 90°: distance horizontal, leaves vertical).
    icoord = np.array(dn["icoord"])
    dcoord = np.array(dn["dcoord"])
    for irow, drow in zip(icoord, dcoord):
        y_scaled = (irow - 5) / 10
        ax_dendro.plot(drow, y_scaled, color="#444", linewidth=1.0)
    ax_dendro.set_ylim(n - 0.5, -0.5)
    ax_dendro.set_xlim(dcoord.max() * 1.05, 0)
    ax_dendro.set_xticks([]); ax_dendro.set_yticks([])
    for sp in ax_dendro.spines.values(): sp.set_visible(False)

    # Annotation strip with sample name overlay.
    label_font = max(6, min(11, 1100 // n))
    for i, s in enumerate(samples_ord):
        batch = sample_batch_map.get(s, "")
        col = _hex_for_batch(batch or "n/a", batch_colors) if batch else "#dddddd"
        ax_annot.add_patch(mpatches.Rectangle(
            (-0.5, i - 0.5), 1, 1, facecolor=col, edgecolor="none"))
        lum = _hex_to_lum(col)
        txt_col = "white" if lum < 0.55 else "#1a1a1a"
        ax_annot.text(0, i, s, fontsize=label_font - 1,
                      ha="center", va="center", color=txt_col)
    ax_annot.set_xlim(-0.5, 0.5); ax_annot.set_ylim(n - 0.5, -0.5)
    ax_annot.set_xticks([]); ax_annot.set_yticks([])
    for sp in ax_annot.spines.values(): sp.set_visible(False)

    # Heatmap lower triangle (with optional in-cell value labels for small n).
    cell_font = max(5, min(9, 700 // n))
    do_labels = n <= 50
    for i in range(n):
        for j in range(n):
            if i <= j: continue
            v = int(M_ord[i, j])
            norm_v = norm(min(v, vmax))
            color = SNP_CMAP(norm_v)
            ax_heat.add_patch(mpatches.Rectangle(
                (j - 0.5, i - 0.5), 1, 1, facecolor=color,
                edgecolor="#cccccc", linewidth=0.2))
            if do_labels:
                txt = "white" if norm_v > 0.58 else "#1a1a1a"
                ax_heat.text(j, i, str(v), ha="center", va="center",
                             fontsize=cell_font, color=txt)
    ax_heat.set_xlim(-0.5, n - 0.5); ax_heat.set_ylim(n - 0.5, -0.5)
    ax_heat.set_xticks(range(n))
    # X labels at 45° (per user request).
    ax_heat.set_xticklabels(samples_ord, rotation=45, ha="right",
                            fontsize=label_font)
    # No Y labels on the heatmap itself — sample IDs are already in the strip
    # on the left (annotation column overlaid with sample name + batch colour).
    ax_heat.set_yticks([])
    for sp in ax_heat.spines.values(): sp.set_visible(False)

    cax = fig.add_axes([0.88, 0.20, 0.012, 0.35])
    sm = plt.cm.ScalarMappable(cmap=SNP_CMAP, norm=norm)
    cb = plt.colorbar(sm, cax=cax)
    cb.set_label(f"Pair-wise SNP distance (capped at {vmax})",
                 fontsize=label_font + 1)
    cb.ax.tick_params(labelsize=label_font)

    fig.suptitle(f"SNP distance heatmap — {cohort} (n={n}) · UPGMA average linkage",
                 fontsize=label_font + 4, weight="bold", y=0.97)
    return f'<img class="fig" src="{_fig_to_b64(fig, dpi=140)}" alt="snp_heat_{cohort}">'


def _render_msn_upgma(samples: list[str], M: np.ndarray, cohort: str,
                      sample_batch_map: dict, batch_colors: dict,
                      upgma_threshold: float = 12.0,
                      export_dir: Optional[Path] = None) -> str:
    """Minimum Spanning Network with UPGMA cluster halos.
    Node colour ← report's batch palette (consistent with the rest of the
    document).  Layout via graphviz neato if available, fallback to spring.

    When `export_dir` is given, the figure is also saved as a standalone
    high-resolution PNG (220 DPI) for direct download.  A link to it is
    appended to the returned HTML.
    """
    import shutil, subprocess, tempfile
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import minimum_spanning_tree
    from scipy.spatial import ConvexHull
    import matplotlib.patches as mpatches
    import networkx as nx

    n = len(samples)
    if n < 3: return ""

    mat = M.astype(float).copy()
    np.fill_diagonal(mat, np.inf)
    mst_arr = minimum_spanning_tree(csr_matrix(mat)).toarray()
    G = nx.Graph(); G.add_nodes_from(samples)
    mst_pairs = set()
    for i in range(n):
        for j in range(n):
            if mst_arr[i, j] > 0:
                a, b = samples[i], samples[j]
                G.add_edge(a, b, snp=int(mst_arr[i, j]), kind="mst")
                mst_pairs.add(tuple(sorted([a, b])))
    # Tie edges: alternative shortest paths with the same SNP weight as the
    # MST edge incident on either endpoint.  Adds the dashed light-blue links
    # that show "this pair would also be equally close in another topology".
    for (a, b) in list(mst_pairs):
        w = G[a][b]["snp"]
        for node in samples:
            if node in (a, b): continue
            for ep in (a, b):
                ep_idx = samples.index(ep); n_idx = samples.index(node)
                if int(M[ep_idx, n_idx]) == w:
                    pair = tuple(sorted([ep, node]))
                    if pair not in mst_pairs and not G.has_edge(*pair):
                        G.add_edge(ep, node, snp=int(w), kind="tie")

    # Layout: prefer graphviz neato for large graphs.  Larger sep="+25" and
    # bigger edge lengths so node labels and SNP labels stay readable.
    pos = None
    if shutil.which("neato"):
        try:
            lines = ["graph G {", "  overlap=false;", '  splines=true;', '  sep="+40";']
            for nm in G.nodes():
                lines.append(f'  "{nm}";')
            for u, v, d in G.edges(data=True):
                length = 0.8 + min(d["snp"], 80) * 0.06
                lines.append(f'  "{u}" -- "{v}" [len={length:.3f}];')
            lines.append("}")
            with tempfile.NamedTemporaryFile("w", suffix=".dot", delete=False) as f:
                f.write("\n".join(lines)); dot_path = f.name
            res = subprocess.run(["neato", "-Tplain", dot_path],
                                 capture_output=True, text=True, timeout=120)
            if res.returncode == 0:
                pos = {}
                for line in res.stdout.splitlines():
                    parts = line.split()
                    if parts and parts[0] == "node":
                        pos[parts[1].strip('"')] = (float(parts[2]), float(parts[3]))
            Path(dot_path).unlink(missing_ok=True)
        except Exception:
            pos = None
    if not pos:
        pos = nx.spring_layout(G, seed=42, k=0.6, iterations=400)

    # UPGMA clusters + halos.
    Z = linkage(squareform(M.astype(float), checks=False), method="average")
    clusters = dict(zip(samples, fcluster(Z, t=upgma_threshold, criterion="distance")))
    clonality = {}
    for i, s in enumerate(samples):
        row = np.delete(M[i].astype(float), i)
        nn = np.sort(row)[:5] if len(row) else np.array([0.0])
        clonality[s] = 1.0 / (float(np.median(nn)) + 1.0)

    # Normalise coordinates to [0, 1].
    xs = np.array([p[0] for p in pos.values()])
    ys = np.array([p[1] for p in pos.values()])
    rx, ry = max(xs.max() - xs.min(), 1e-9), max(ys.max() - ys.min(), 1e-9)
    pos_n = {k: ((v[0] - xs.min()) / rx, (v[1] - ys.min()) / ry)
             for k, v in pos.items()}

    # Larger canvas + higher DPI on render so node labels and SNP edge values
    # are clearly readable even with ~60 nodes.  Extra axis padding so cluster
    # halos near the edges don't get clipped at the image border.
    fig = plt.figure(figsize=(18, 13))
    ax = fig.add_axes([0.02, 0.05, 0.72, 0.92])
    ax.axis("off")
    ax.set_xlim(-0.18, 1.18); ax.set_ylim(-0.18, 1.18); ax.set_aspect("equal")

    # Cluster halos (convex hull) for groups of ≥3 samples.
    HALO_PALETTE = ["#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5",
                    "#bc80bd", "#ccebc5", "#ffed6f", "#8dd3c7", "#bebada"]
    by_cluster: dict = {}
    for nm, cid in clusters.items():
        by_cluster.setdefault(cid, []).append(nm)
    big = {cid: ms for cid, ms in by_cluster.items() if len(ms) >= 3}
    halo_colour = {cid: HALO_PALETTE[i % len(HALO_PALETTE)]
                   for i, cid in enumerate(sorted(big))}
    for cid, members in big.items():
        pts = np.array([pos_n[m] for m in members])
        if len(pts) >= 3:
            try:
                hull = ConvexHull(pts)
                hp = pts[hull.vertices]
                cent = hp.mean(axis=0)
                # Expansion factor reduced 1.40 → 1.25 so halos stay closer to
                # their cluster members and don't push past the axis padding.
                exp = cent + (hp - cent) * 1.25
                ax.fill(exp[:, 0], exp[:, 1], color=halo_colour[cid], alpha=0.30,
                        zorder=0, edgecolor=halo_colour[cid], linewidth=2.0,
                        linestyle="--")
            except Exception:
                pass

    # Edges — thin solid grey for MST, dashed light-blue for ties.
    for u, v, d in G.edges(data=True):
        x1, y1 = pos_n[u]; x2, y2 = pos_n[v]
        snp = d["snp"]
        if d.get("kind") == "tie":
            colour, ls, lw = "#5B9ABF", "--", 0.9
        else:
            colour, ls, lw = "#555", "-", 0.7
        ax.plot([x1, x2], [y1, y2], color=colour, linewidth=lw, linestyle=ls,
                alpha=0.7, zorder=1)
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        if snp <= 5:    fs, weight = 9, "bold"
        elif snp <= 15: fs, weight = 8, "bold"
        else:           fs, weight = 7, "normal"
        ax.text(mx, my, str(snp), fontsize=fs, weight=weight,
                ha="center", va="center", zorder=3,
                bbox=dict(boxstyle="round,pad=0.12", facecolor="white",
                          edgecolor="#ccc", linewidth=0.3, alpha=0.92))

    # Nodes — constant size, colour from batch palette (consistent with report).
    NODE_SIZE = 380
    for nm in G.nodes():
        batch = sample_batch_map.get(nm, "")
        col = _hex_for_batch(batch or "n/a", batch_colors) if batch else "#bbbbbb"
        ax.scatter(*pos_n[nm], s=NODE_SIZE, c=[col],
                   edgecolors="#222", linewidths=1.0, zorder=4)
        ax.text(pos_n[nm][0], pos_n[nm][1] - 0.020, nm, fontsize=7,
                ha="center", va="top", weight="600", zorder=5,
                bbox=dict(boxstyle="round,pad=0.12", facecolor="white",
                          edgecolor="#aaa", linewidth=0.3, alpha=0.92))

    title = (f"MSN + UPGMA — {cohort} (n={n} · cut @ {upgma_threshold:g} SNPs"
             f" → {len(big)} clusters)")
    ax.text(0.5, 1.01, title, transform=ax.transAxes, fontsize=13,
            weight="bold", ha="center", va="bottom")

    # Side legends — batch only if there is more than one or any batch labelled.
    batches_seen = sorted({sample_batch_map.get(s, "") for s in samples
                           if sample_batch_map.get(s, "")})
    if batches_seen:
        ax_leg = fig.add_axes([0.77, 0.50, 0.22, 0.42])
        ax_leg.axis("off")
        ax_leg.set_title("Batch", fontsize=12, weight="bold", loc="left", pad=4)
        rh = 0.92 / max(len(batches_seen), 1)
        for i, b in enumerate(batches_seen):
            y = 0.93 - i * rh
            col = _hex_for_batch(b, batch_colors)
            cnt = sum(1 for s in samples if sample_batch_map.get(s) == b)
            ax_leg.scatter([0.06], [y], s=240, c=[col],
                           edgecolors="#222", linewidths=1.2,
                           transform=ax_leg.transAxes)
            ax_leg.text(0.20, y, f"{b}  (n={cnt})", fontsize=10, va="center",
                        transform=ax_leg.transAxes)
    if big:
        ax_cl = fig.add_axes([0.77, 0.05, 0.22, 0.40])
        ax_cl.axis("off")
        ax_cl.set_title(f"UPGMA clusters @ {upgma_threshold:g} SNPs",
                        fontsize=12, weight="bold", loc="left", pad=4)
        top = sorted(big.items(), key=lambda kv: -len(kv[1]))[:8]
        rh = 0.92 / max(len(top), 1)
        for i, (cid, members) in enumerate(top):
            y = 0.93 - i * rh
            col = halo_colour[cid]
            ax_cl.add_patch(mpatches.Rectangle(
                (0.04, y - 0.022), 0.10, 0.04, facecolor=col, alpha=0.5,
                edgecolor=col, linewidth=1.6, linestyle="--",
                transform=ax_cl.transAxes))
            ax_cl.text(0.18, y, f"C{cid} · n={len(members)}", fontsize=10,
                       va="center", transform=ax_cl.transAxes)

    # Save a standalone high-DPI PNG sibling to the HTML for direct download.
    export_link = ""
    if export_dir is not None:
        try:
            export_dir.mkdir(parents=True, exist_ok=True)
            png_path = export_dir / f"msn_upgma_{cohort}.png"
            fig.savefig(png_path, dpi=220, bbox_inches="tight", facecolor="white")
            export_link = (
                f'<p style="margin:6px 0 4px;font-size:11px;">'
                f'<a href="{png_path.name}" target="_blank" '
                f'style="color:{BLUE};font-weight:600;">'
                f'Download high-resolution PNG ({png_path.name}) →</a></p>')
        except Exception:
            pass
    return export_link + f'<img class="fig" src="{_fig_to_b64(fig, dpi=180)}" alt="msn_{cohort}">'


_SSN_COUNTER = 0


def _render_ssn_slider(samples: list[str], M: np.ndarray, cohort: str,
                       sample_batch_map: dict, batch_colors: dict) -> str:
    """Interactive Sample Similarity Network with D3.js force layout.

    Behaviour modelled on the make_network.py template:
      - Drag any node to move it; on release the force simulation pulls it
        back to its energy-minimal position (the "spring-back" effect).
      - SNP threshold slider rebuilds the edge set live.
      - Click a node to see its metadata and highlight its incident edges.
      - Click the batch legend to filter samples by hospital / batch.
      - Hover over a node for tooltip with id + batch + collection_date.
      - Zoom & pan with mouse-wheel / drag on empty canvas.

    Node colours come from the same batch palette used elsewhere in the
    report so HGRAL/HLAFE (or whichever batches are present) stay consistent.
    """
    global _SSN_COUNTER
    n = len(samples)
    if n < 3:
        return ""
    _SSN_COUNTER += 1
    uid = f"ssn_{_SSN_COUNTER}"

    # Pre-compute every pairwise edge with its SNP distance — D3 filters in JS.
    edges_all = []
    for i in range(n):
        for j in range(i + 1, n):
            edges_all.append({"source": samples[i],
                              "target": samples[j],
                              "snp": int(M[i, j])})
    max_snp = max((e["snp"] for e in edges_all), default=0)
    # Default threshold = 30th percentile (informative but not saturated).
    p30 = int(np.percentile([e["snp"] for e in edges_all], 30)) if edges_all else 0
    initial_thr = max(1, p30)

    nodes = []
    for s in samples:
        b = sample_batch_map.get(s, "") or "Unknown"
        nodes.append({
            "id": s,
            "batch": b,
            "color": _hex_for_batch(b if b != "Unknown" else "n/a", batch_colors)
                     if b != "Unknown" else "#bbbbbb",
        })

    data_json = json.dumps({"nodes": nodes, "edges": edges_all})

    # The CSS is scoped to `#{uid}` so multiple SSNs can coexist in one page.
    return (
        f'<div id="{uid}" class="ssn-wrap">'
        f'  <style>'
        f'    #{uid} {{ display:flex; height:720px; border:1px solid #e5e7eb; '
        f'              border-radius:8px; overflow:hidden; background:#fff; }}'
        f'    #{uid} .ssn-canvas {{ flex:1; position:relative; background:#fff; }}'
        f'    #{uid} .ssn-side   {{ width:260px; border-left:1px solid #e5e7eb; '
        f'                          padding:14px; overflow-y:auto; font-size:12px; '
        f'                          background:#FAFBFC; }}'
        f'    #{uid} .ssn-side h4{{ margin:10px 0 6px; font-size:11px; color:#555; '
        f'                          text-transform:uppercase; letter-spacing:0.04em; '
        f'                          font-weight:700; }}'
        f'    #{uid} .ssn-slider-box {{ background:#fff; padding:10px; '
        f'                              border-radius:6px; border:1px solid #e5e7eb; }}'
        f'    #{uid} .ssn-slider-box input[type=range] {{ width:100%; accent-color:{BLUE}; }}'
        f'    #{uid} .ssn-slider-box .v {{ font-size:18px; font-weight:700; '
        f'                                  color:{ORANGE}; text-align:center; }}'
        f'    #{uid} .ssn-stat {{ margin-bottom:3px; }}'
        f'    #{uid} .ssn-stat b {{ color:#222; }}'
        f'    #{uid} .ssn-leg-item {{ display:flex; align-items:center; padding:3px 2px; '
        f'                            cursor:pointer; border-radius:3px; }}'
        f'    #{uid} .ssn-leg-item:hover {{ background:#eef; }}'
        f'    #{uid} .ssn-leg-item.off {{ opacity:0.3; }}'
        f'    #{uid} .ssn-leg-dot {{ width:13px; height:13px; border-radius:50%; '
        f'                           margin-right:8px; border:1.5px solid #333; }}'
        f'    #{uid} .ssn-leg-cnt {{ margin-left:auto; color:#888; font-size:11px; }}'
        f'    #{uid} .ssn-node {{ stroke:#222; stroke-width:1.2px; cursor:pointer; }}'
        f'    #{uid} .ssn-node:hover {{ stroke-width:3px; }}'
        f'    #{uid} .ssn-node.selected {{ stroke:{ORANGE}; stroke-width:4px; }}'
        f'    #{uid} .ssn-node.dimmed {{ opacity:0.2; }}'
        f'    #{uid} .ssn-link {{ stroke-opacity:0.7; }}'
        f'    #{uid} .ssn-link.highlight {{ stroke:{ORANGE} !important; '
        f'                                   stroke-opacity:1; stroke-width:2.2px !important; }}'
        f'    #{uid} .ssn-label {{ font-size:9px; fill:#333; pointer-events:none; '
        f'                          text-anchor:middle; }}'
        f'    #{uid} .ssn-elabel {{ font-size:8px; pointer-events:none; '
        f'                           text-anchor:middle; font-weight:700; '
        f'                           paint-order:stroke; stroke:#fff; stroke-width:2.5px; '
        f'                           stroke-linejoin:round; }}'
        f'    #{uid} .ssn-snp-legend {{ display:flex; flex-direction:column; gap:3px; '
        f'                              font-size:11px; margin-top:4px; }}'
        f'    #{uid} .ssn-snp-bin {{ display:flex; align-items:center; padding:3px 4px; '
        f'                            cursor:pointer; border-radius:3px; }}'
        f'    #{uid} .ssn-snp-bin:hover {{ background:#eef; }}'
        f'    #{uid} .ssn-snp-bin.off {{ opacity:0.3; }}'
        f'    #{uid} .ssn-snp-sw  {{ width:18px; height:4px; border-radius:2px; margin-right:6px; }}'
        f'    #{uid} .ssn-snp-cnt {{ margin-left:auto; color:#888; font-size:10px; }}'
        f'    #{uid} .ssn-tip {{ position:absolute; padding:6px 8px; '
        f'                       background:rgba(33,37,41,0.94); color:#fff; '
        f'                       border-radius:4px; font-size:11px; pointer-events:none; '
        f'                       opacity:0; transition:opacity 120ms; z-index:50; }}'
        f'    #{uid} .ssn-info {{ background:#fff; border:1px solid #e5e7eb; '
        f'                        padding:8px; border-radius:6px; font-size:11px; }}'
        f'    #{uid} .ssn-btn {{ background:{BLUE}; color:#fff; border:none; '
        f'                       padding:4px 10px; border-radius:14px; cursor:pointer; '
        f'                       font-size:11px; margin-right:4px; }}'
        f'    #{uid} .ssn-btn:hover {{ background:#1d4364; }}'
        f'  </style>'
        f'  <div class="ssn-canvas">'
        f'    <svg id="{uid}_svg" width="100%" height="100%"></svg>'
        f'    <div id="{uid}_tip" class="ssn-tip"></div>'
        f'  </div>'
        f'  <div class="ssn-side">'
        f'    <h4>SNP threshold</h4>'
        f'    <div class="ssn-slider-box">'
        f'      <div>edges with SNP ≤ <b id="{uid}_thr_val" class="v">{initial_thr}</b></div>'
        f'      <input id="{uid}_thr" type="range" min="0" max="{max_snp}" value="{initial_thr}">'
        f'      <div style="display:flex;justify-content:space-between;font-size:10px;color:#888;">'
        f'        <span>0</span><span>{max_snp}</span></div>'
        f'    </div>'
        f'    <h4>Stats</h4>'
        f'    <div class="ssn-stat">Nodes shown: <b id="{uid}_s_n">—</b></div>'
        f'    <div class="ssn-stat">Edges shown: <b id="{uid}_s_e">—</b></div>'
        f'    <div class="ssn-stat">Components:  <b id="{uid}_s_c">—</b></div>'
        f'    <div class="ssn-stat">Largest comp: <b id="{uid}_s_lc">—</b></div>'
        f'    <h4>Node repulsion</h4>'
        f'    <div class="ssn-slider-box">'
        f'      <div>charge: <b id="{uid}_spr_val">100</b></div>'
        f'      <input id="{uid}_spread" type="range" min="20" max="400" step="10" value="100">'
        f'      <div style="font-size:10px;color:#888;margin-top:2px;">'
        f'        Lower = nodes pulled tighter into the canvas.</div>'
        f'    </div>'
        f'    <h4>SNP edge scale</h4>'
        f'    <div class="ssn-slider-box">'
        f'      <div>length factor: <b id="{uid}_lks_val">18</b></div>'
        f'      <input id="{uid}_lks" type="range" min="5" max="80" step="1" value="18">'
        f'      <div style="font-size:10px;color:#888;margin-top:2px;">'
        f'        Higher = edges grow more with their SNP count (so distant pairs '
        f'        sit further apart). Edge length = log(SNPs+1) &times; factor.</div>'
        f'    </div>'
        f'    <h4>Controls</h4>'
        f'    <button class="ssn-btn" id="{uid}_btn_labels">sample labels</button>'
        f'    <button class="ssn-btn" id="{uid}_btn_edge">edge labels</button>'
        f'    <button class="ssn-btn" id="{uid}_btn_focus">focus selected</button>'
        f'    <button class="ssn-btn" id="{uid}_btn_reset">reset zoom</button>'
        f'    <h4>Edge colour (SNPs) — click to toggle</h4>'
        f'    <div class="ssn-snp-legend" id="{uid}_bin_legend">'
        f'      <div class="ssn-snp-bin"  data-bin="b1" data-color="#000000">'
        f'        <span class="ssn-snp-sw" style="background:#000000"></span>'
        f'        0–1 SNPs<span class="ssn-snp-cnt" id="{uid}_cnt_b1">—</span></div>'
        f'      <div class="ssn-snp-bin"  data-bin="b2" data-color="#c92020">'
        f'        <span class="ssn-snp-sw" style="background:#c92020"></span>'
        f'        2–3<span class="ssn-snp-cnt" id="{uid}_cnt_b2">—</span></div>'
        f'      <div class="ssn-snp-bin"  data-bin="b3" data-color="#2ca02c">'
        f'        <span class="ssn-snp-sw" style="background:#2ca02c"></span>'
        f'        4–11<span class="ssn-snp-cnt" id="{uid}_cnt_b3">—</span></div>'
        f'      <div class="ssn-snp-bin"  data-bin="b4" data-color="#5B9ABF">'
        f'        <span class="ssn-snp-sw" style="background:#5B9ABF"></span>'
        f'        &ge; 12<span class="ssn-snp-cnt" id="{uid}_cnt_b4">—</span></div>'
        f'    </div>'
        f'    <h4>Batch (click to filter)</h4>'
        f'    <div id="{uid}_legend"></div>'
        f'    <h4>Selected sample</h4>'
        f'    <div class="ssn-info" id="{uid}_info">click a node…</div>'
        f'  </div>'
        f'</div>'
        f'<script>'
        f'(function(){{'
        f'  const DATA = {data_json};'
        f'  const UID = "{uid}";'
        f'  const wrap = document.getElementById("{uid}");'
        f'  const canvas = wrap.querySelector(".ssn-canvas");'
        f'  const W = canvas.clientWidth || 800, H = canvas.clientHeight || 720;'
        f'  const svg = d3.select("#{uid}_svg").attr("viewBox",[0,0,W,H]);'
        f'  const g   = svg.append("g");'
        f'  const zoom = d3.zoom().scaleExtent([0.1,8]).on("zoom",e=>g.attr("transform",e.transform));'
        f'  svg.call(zoom);'
        f'  const linkG  = g.append("g").attr("class","links");'
        f'  const elblG  = g.append("g").attr("class","elabels");'
        f'  const nodeG  = g.append("g").attr("class","nodes");'
        f'  const lblG   = g.append("g").attr("class","labels");'
        f'  const nMap  = new Map(DATA.nodes.map(d => [d.id, {{...d}}]));'
        f'  const nodes = [...nMap.values()];'
        f'  nodes.forEach(d => {{ d.x = W/2 + (Math.random()-0.5)*200; d.y = H/2 + (Math.random()-0.5)*200; }});'
        f'  const tip = d3.select("#{uid}_tip");'
        f'  let active = []; let curThr = {initial_thr};'
        f'  let enabled = new Set(DATA.nodes.map(d => d.batch));'
        f'  let labelsOn = true;'
        f'  let edgeLabelsOn = true;'
        f'  let charge = -100;'
        f'  let snpScale = 18;'
        f'  /* link distance grows with SNP count via log scaling — pairs with'
        f'     many SNPs sit further apart, clonal pairs stay tight. */'
        f'  const linkDist = d => Math.max(15, Math.log(d.snp + 1) * snpScale);'
        f'  /* Edge colour by SNP bin (0-1 black · 2-3 red · 4-11 green · ≥12 blue). */'
        f'  function snpBin(snp){{'
        f'    if (snp <=  1) return "b1";'
        f'    if (snp <=  3) return "b2";'
        f'    if (snp <= 11) return "b3";'
        f'    return "b4";'
        f'  }}'
        f'  const BIN_COLOUR = {{ b1:"#000000", b2:"#c92020", b3:"#2ca02c", b4:"#5B9ABF" }};'
        f'  function snpColour(snp){{ return BIN_COLOUR[snpBin(snp)]; }}'
        f'  /* Per-bin show/hide.  All bins ON by default — at threshold=max the graph'
        f'     must show every pair.  User can click any legend row to declutter. */'
        f'  let binsOn = {{ b1:true, b2:true, b3:true, b4:true }};'
        f'  /* Focus mode: when ON and a node is selected, only its incident'
        f'     edges (and itself) are drawn — everything else is hidden. */'
        f'  let focusOn = false; let focusNode = null;'
        f'  /* forceX / forceY constrain the cloud to the canvas; charge controls'
        f'     how strongly nodes repel.  Bounded layout stops nodes flying off. */'
        f'  const sim = d3.forceSimulation(nodes)'
        f'    .force("link", d3.forceLink([]).id(d=>d.id).distance(linkDist).strength(0.7))'
        f'    .force("charge", d3.forceManyBody().strength(-100))'
        f'    .force("center", d3.forceCenter(W/2, H/2))'
        f'    .force("x", d3.forceX(W/2).strength(0.07))'
        f'    .force("y", d3.forceY(H/2).strength(0.07))'
        f'    .force("collide", d3.forceCollide(14));'
        f'  function updateNet(thr){{'
        f'    /* 1) Filter by SNP threshold + bin toggles + batch toggles + focus. */'
        f'    active = DATA.edges.filter(e => e.snp <= thr'
        f'      && binsOn[snpBin(e.snp)]'
        f'      && enabled.has(nMap.get(e.source).batch)'
        f'      && enabled.has(nMap.get(e.target).batch)'
        f'      && (!focusOn || !focusNode'
        f'           || e.source === focusNode || e.target === focusNode'
        f'           || (typeof e.source === "object" && e.source.id === focusNode)'
        f'           || (typeof e.target === "object" && e.target.id === focusNode)))'
        f'      .map(e => ({{...e}}));'
        f'    /* 2) Per-bin counters for the legend. */'
        f'    const binCnt = {{ b1:0, b2:0, b3:0, b4:0 }};'
        f'    DATA.edges.forEach(e => {{ if (e.snp <= thr) binCnt[snpBin(e.snp)]++; }});'
        f'    document.getElementById(UID+"_cnt_b1").textContent = binCnt.b1;'
        f'    document.getElementById(UID+"_cnt_b2").textContent = binCnt.b2;'
        f'    document.getElementById(UID+"_cnt_b3").textContent = binCnt.b3;'
        f'    document.getElementById(UID+"_cnt_b4").textContent = binCnt.b4;'
        f'    const link = linkG.selectAll("line").data(active, d => d.source+"-"+d.target);'
        f'    link.exit().remove();'
        f'    const linkEnter = link.enter().append("line").attr("class","ssn-link")'
        f'      .attr("stroke-width", 0.7)'
        f'      .attr("stroke", d => snpColour(d.snp));'
        f'    linkEnter.append("title").text(d => d.snp + " SNPs");'
        f'    /* Re-paint existing edges too in case the threshold change re-uses '
        f'       previously-rendered <line>s with a different SNP value. */'
        f'    linkG.selectAll("line")'
        f'      .attr("stroke", d => snpColour(d.snp))'
        f'      .attr("stroke-width", 0.7);'
        f'    /* Edge labels: append/remove SVG text for each visible edge. */'
        f'    const el = elblG.selectAll("text").data(active, d => d.source+"-"+d.target);'
        f'    el.exit().remove();'
        f'    el.enter().append("text").attr("class","ssn-elabel")'
        f'      .attr("fill", d => snpColour(d.snp))'
        f'      .text(d => d.snp);'
        f'    elblG.selectAll("text")'
        f'      .attr("fill", d => snpColour(d.snp))'
        f'      .text(d => d.snp);'
        f'    elblG.style("display", edgeLabelsOn ? null : "none");'
        f'    sim.force("link").links(active);'
        f'    sim.alpha(0.4).restart();'
        f'    const vNodes = nodes.filter(d => enabled.has(d.batch));'
        f'    document.getElementById(UID+"_s_n").textContent = vNodes.length;'
        f'    document.getElementById(UID+"_s_e").textContent = active.length;'
        f'    const adj = new Map(vNodes.map(n=>[n.id,new Set()]));'
        f'    active.forEach(e=>{{'
        f'      const sId = e.source.id||e.source, tId = e.target.id||e.target;'
        f'      if (adj.has(sId)) adj.get(sId).add(tId);'
        f'      if (adj.has(tId)) adj.get(tId).add(sId);'
        f'    }});'
        f'    const vis = new Set(); const sizes = [];'
        f'    for (const v of vNodes){{ if (vis.has(v.id)) continue;'
        f'      let q=[v.id], sz=0; while(q.length){{ const c=q.pop();'
        f'        if (vis.has(c)) continue; vis.add(c); sz++;'
        f'        for (const nb of (adj.get(c)||[])) if (!vis.has(nb)) q.push(nb); }}'
        f'      sizes.push(sz); }}'
        f'    document.getElementById(UID+"_s_c").textContent  = sizes.length;'
        f'    document.getElementById(UID+"_s_lc").textContent = sizes.length ? Math.max(...sizes) : 0;'
        f'  }}'
        f'  /* Render nodes once; drag uses spring-back via release of fx/fy. */'
        f'  const node = nodeG.selectAll("circle").data(nodes, d=>d.id).enter().append("circle")'
        f'    .attr("class","ssn-node").attr("r",7).attr("fill",d=>d.color);'
        f'  node.call(d3.drag()'
        f'    .on("start",(e,d)=>{{ if (!e.active) sim.alphaTarget(0.3).restart(); d.fx=d.x; d.fy=d.y; }})'
        f'    .on("drag",(e,d)=>{{ d.fx=e.x; d.fy=e.y; }})'
        f'    .on("end",(e,d)=>{{ if (!e.active) sim.alphaTarget(0); d.fx=null; d.fy=null; }}));'
        f'  node.on("mouseover",(e,d)=>{{'
        f'    tip.style("opacity",1).html("<b>"+d.id+"</b><br>batch: "+d.batch); }})'
        f'   .on("mousemove",(e)=>{{ tip.style("left",(e.offsetX+14)+"px").style("top",(e.offsetY+14)+"px"); }})'
        f'   .on("mouseout",()=>tip.style("opacity",0))'
        f'   .on("click",(e,d)=>{{'
        f'     nodeG.selectAll("circle").classed("selected",false);'
        f'     d3.select(e.currentTarget).classed("selected",true);'
        f'     focusNode = d.id;'
        f'     document.getElementById(UID+"_info").innerHTML ='
        f'       "<b>"+d.id+"</b><br>batch: "+d.batch;'
        f'     linkG.selectAll("line").classed("highlight",'
        f'       l => (l.source.id||l.source)===d.id || (l.target.id||l.target)===d.id);'
        f'     if (focusOn) updateNet(curThr);'
        f'   }});'
        f'  const lbl = lblG.selectAll("text").data(nodes, d=>d.id).enter().append("text")'
        f'    .attr("class","ssn-label").attr("dy",-10).text(d=>d.id);'
        f'  sim.on("tick", ()=>{{'
        f'    /* Clamp positions inside the canvas so nodes never fly out. */'
        f'    const PAD = 25;'
        f'    nodes.forEach(d => {{ d.x = Math.max(PAD, Math.min(W-PAD, d.x));'
        f'                          d.y = Math.max(PAD, Math.min(H-PAD, d.y)); }});'
        f'    linkG.selectAll("line")'
        f'      .attr("x1",d=>d.source.x).attr("y1",d=>d.source.y)'
        f'      .attr("x2",d=>d.target.x).attr("y2",d=>d.target.y);'
        f'    elblG.selectAll("text")'
        f'      .attr("x",d=>(d.source.x+d.target.x)/2)'
        f'      .attr("y",d=>(d.source.y+d.target.y)/2 - 2);'
        f'    nodeG.selectAll("circle").attr("cx",d=>d.x).attr("cy",d=>d.y);'
        f'    lblG.selectAll("text").attr("x",d=>d.x).attr("y",d=>d.y);'
        f'  }});'
        f'  /* Slider + buttons + batch legend. */'
        f'  const sld = document.getElementById(UID+"_thr");'
        f'  const sldVal = document.getElementById(UID+"_thr_val");'
        f'  sld.addEventListener("input", e=>{{ curThr=+e.target.value; sldVal.textContent=curThr; updateNet(curThr); }});'
        f'  document.getElementById(UID+"_btn_labels").onclick = ()=>{{ labelsOn=!labelsOn; lblG.style("display",labelsOn?null:"none"); }};'
        f'  document.getElementById(UID+"_btn_edge").onclick   = ()=>{{ edgeLabelsOn=!edgeLabelsOn; elblG.style("display",edgeLabelsOn?null:"none"); }};'
        f'  document.getElementById(UID+"_btn_focus").onclick  = function(){{'
        f'    focusOn = !focusOn;'
        f'    this.style.background = focusOn ? "{ORANGE}" : "{BLUE}";'
        f'    updateNet(curThr);'
        f'  }};'
        f'  document.getElementById(UID+"_btn_reset").onclick  = ()=>{{ svg.transition().duration(400).call(zoom.transform, d3.zoomIdentity); }};'
        f'  /* Bin legend click-to-toggle: each row in the SNP colour legend can'
        f'     hide its edges so the user can declutter dense graphs. */'
        f'  document.querySelectorAll("#"+UID+"_bin_legend .ssn-snp-bin").forEach(el => {{'
        f'    el.addEventListener("click", () => {{'
        f'      const k = el.getAttribute("data-bin");'
        f'      binsOn[k] = !binsOn[k];'
        f'      el.classList.toggle("off", !binsOn[k]);'
        f'      updateNet(curThr);'
        f'    }});'
        f'  }});'
        f'  /* Node-repulsion slider: tune forceManyBody charge. */'
        f'  const spr = document.getElementById(UID+"_spread");'
        f'  const sprVal = document.getElementById(UID+"_spr_val");'
        f'  spr.addEventListener("input", e=>{{'
        f'    charge = -(+e.target.value);'
        f'    sprVal.textContent = (+e.target.value);'
        f'    sim.force("charge").strength(charge);'
        f'    sim.alpha(0.4).restart();'
        f'  }});'
        f'  /* SNP edge-scale slider: tune how much SNP count stretches each edge. */'
        f'  const lks = document.getElementById(UID+"_lks");'
        f'  const lksVal = document.getElementById(UID+"_lks_val");'
        f'  lks.addEventListener("input", e=>{{'
        f'    snpScale = +e.target.value;'
        f'    lksVal.textContent = snpScale;'
        f'    sim.force("link").distance(linkDist);'
        f'    sim.alpha(0.5).restart();'
        f'  }});'
        f'  const cnt = {{}}; DATA.nodes.forEach(n => cnt[n.batch] = (cnt[n.batch]||0)+1);'
        f'  const leg = d3.select("#"+UID+"_legend");'
        f'  Object.entries(cnt).sort((a,b)=>b[1]-a[1]).forEach(([b,c])=>{{'
        f'    const it = leg.append("div").attr("class","ssn-leg-item").on("click", function(){{'
        f'      if (enabled.has(b)) {{ enabled.delete(b); d3.select(this).classed("off",true); }}'
        f'      else                 {{ enabled.add(b);    d3.select(this).classed("off",false); }}'
        f'      nodeG.selectAll("circle").classed("dimmed", d => !enabled.has(d.batch));'
        f'      lblG.selectAll("text").classed("dimmed",  d => !enabled.has(d.batch));'
        f'      updateNet(curThr);'
        f'    }});'
        f'    const colour = DATA.nodes.find(n => n.batch === b).color;'
        f'    it.append("div").attr("class","ssn-leg-dot").style("background", colour);'
        f'    it.append("span").text(b);'
        f'    it.append("span").attr("class","ssn-leg-cnt").text(c);'
        f'  }});'
        f'  svg.on("click",(e)=>{{ if (e.target.tagName==="svg"){{'
        f'    nodeG.selectAll("circle").classed("selected",false);'
        f'    linkG.selectAll("line").classed("highlight",false);'
        f'    focusNode = null;'
        f'    document.getElementById(UID+"_info").innerHTML = "click a node…";'
        f'    if (focusOn) updateNet(curThr);'
        f'  }} }});'
        f'  updateNet(curThr);'
        f'}})();'
        f'</script>'
    )


def render_snp_section(results_dir: Path, rows: list[SampleRow],
                       batch_colors: dict[str, str],
                       export_dir: Optional[Path] = None) -> str:
    """SNP heatmap (Plotly) + MSN/UPGMA (matplotlib) + interactive SSN slider."""
    cohort_dirs = sorted((results_dir / "07_snp" / "snpdists").glob("*/"))
    if not cohort_dirs:
        return ""
    sample_batch_map = {r.sample_id: r.batch for r in rows}

    bits: list[str] = []
    for cd in cohort_dirs:
        cohort = cd.name
        molten = cd / f"{cohort}.snpdist.molten.tsv"
        if not molten.exists():
            continue
        samples, M = _load_snp_matrix(molten)
        if not samples:
            continue
        bits.append(f"<h4>{cohort} — SNP distance heatmap</h4>")
        bits.append('<p class="muted">Lower-triangle pair-wise SNP distance. '
                    'Samples ordered by UPGMA average linkage (left dendrogram); '
                    'batch annotation strip with sample IDs overlaid.</p>')
        bits.append(_render_snp_heatmap_template(samples, M, cohort,
                                                 sample_batch_map, batch_colors))
        bits.append(f"<h4>{cohort} — MSN + UPGMA network</h4>")
        bits.append('<p class="muted">Solid grey = MST edges. '
                    'Dashed light-blue = tied alternative shortest links.</p>')
        bits.append(_render_msn_upgma(samples, M, cohort,
                                      sample_batch_map, batch_colors,
                                      export_dir=export_dir))
        bits.append(f"<h4>{cohort} — SSN (interactive, threshold slider)</h4>")
        bits.append('<p class="muted">Sample Similarity Network with 2-D MDS '
                    'layout (closer in plot ≈ closer in SNPs). Drag the slider '
                    'to widen / tighten the SNP threshold and see which pairs '
                    'connect. Top of range = (almost) complete graph; bottom = '
                    'isolated clusters.</p>')
        bits.append(_render_ssn_slider(samples, M, cohort,
                                       sample_batch_map, batch_colors))
    if not bits:
        return ""
    return f"""
    <section class="section">
      <h3>SNP analysis</h3>
      {''.join(bits)}
    </section>
    """


# ─── Phase 7 — Excel companion ─────────────────────────────────────────────

def write_excel_companion(rows: list[SampleRow],
                          quast_rows: list[dict],
                          busco_rows: list[dict],
                          results_dir: Path,
                          out_xlsx: Path) -> None:
    try:
        import pandas as pd
    except ImportError:
        print("[generate_epicandi_report] pandas missing; skipping Excel.",
              file=sys.stderr)
        return

    out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as xw:

        # Sheet 1 — Summary_batch
        batch_stats = defaultdict(lambda: {"n":0, "pass":0, "amr":0, "cnv":0})
        for r in rows:
            b = r.batch or "(no batch)"
            batch_stats[b]["n"]    += 1
            batch_stats[b]["pass"] += int(r.qc_flag == "PASS")
            batch_stats[b]["amr"]  += int(r.has_amr)
            batch_stats[b]["cnv"]  += int(r.has_cnv)
        pd.DataFrame([{"batch":b, **v} for b,v in sorted(batch_stats.items())])\
          .to_excel(xw, sheet_name="Summary_batch", index=False)

        # Sheet 2 — Summary_species
        sp_stats = Counter(r.species_label for r in rows)
        pd.DataFrame([{"species":k, "samples":v} for k,v in sp_stats.most_common()])\
          .to_excel(xw, sheet_name="Summary_species", index=False)

        # Sheet 3 — Summary_fastq_qc (overall numbers)
        n = len(rows)
        passN = sum(1 for r in rows if r.qc_flag == "PASS")
        clean = sum(r.clean_gb for r in rows)
        pd.DataFrame([{
            "samples_total": n,
            "samples_pass":  passN,
            "samples_fail":  n - passN,
            "clean_GB_total": round(clean, 2),
            "median_est_depth_x":   round(np.median([r.est_depth for r in rows]) if rows else 0, 1),
            "median_align_rate_pct":round(np.median([r.align_rate for r in rows]) if rows else 0, 1),
            "median_sylph_cov_pct": round(np.median([r.sylph_cov for r in rows]) if rows else 0, 1),
        }]).to_excel(xw, sheet_name="Summary_fastq_qc", index=False)

        # Sheet 4 — fastq_qc (per-sample)
        pd.DataFrame([{
            "sample_id":      r.sample_id,
            "platform":       r.platform,
            "batch":          r.batch,
            "is_external":    r.is_external,
            "collection_date":r.collection_date,
            "qc_flag":        r.qc_flag,
            "species":        r.species,
            "clade":          r.clade,
            "reference_slug": r.reference_slug,
            "raw_reads":      r.raw_illu_reads,
            "raw_GB":         round(r.raw_gb, 3),
            "clean_reads":    r.clean_illu_reads,
            "clean_GB":       round(r.clean_gb, 3),
            "retention_pct":  round(r.retention_pct, 2),
            "sylph_ani":      r.sylph_ani,
            "sylph_cov_pct":  r.sylph_cov,
            "align_rate_pct": r.align_rate,
            "est_depth_x":    r.est_depth,
            "contam_pct":     r.contam,
            "qc_notes":       r.qc_notes,
        } for r in rows]).to_excel(xw, sheet_name="fastq_qc", index=False)

        # Sheet 5 — assembly_stats
        if quast_rows:
            pd.DataFrame(quast_rows).to_excel(xw, sheet_name="assembly_stats",
                                              index=False)

        # Sheet 6 — amr_report+cnv (join master_table + per-sample CNV summary)
        master = results_dir / "cohort" / "master_table.tsv"
        if master.exists():
            df_amr = pd.read_csv(master, sep="\t")
        else:
            df_amr = pd.DataFrame()
        cnv_counts = []
        for r in rows:
            cnv_path = results_dir / "08_cnv" / r.sample_id / f"{r.sample_id}.cnv_events.tsv"
            recs = _read_tsv(cnv_path)
            n_evt = sum(1 for x in recs if (x.get("event") or "normal") not in ("normal", ""))
            cnv_counts.append({"sample_id": r.sample_id, "cnv_events": n_evt,
                               "cnv_genes_affected": ",".join(
                                   sorted({x["gene"] for x in recs
                                           if (x.get("event") or "normal") not in ("normal","")})
                               )})
        df_cnv = pd.DataFrame(cnv_counts)
        if not df_amr.empty:
            df_merged = df_amr.merge(df_cnv, left_on=df_amr.columns[0],
                                     right_on="sample_id", how="outer")
        else:
            df_merged = df_cnv
        df_merged.to_excel(xw, sheet_name="amr_report+cnv", index=False)

    print(f"[generate_epicandi_report] wrote {out_xlsx}", file=sys.stderr)


# ═════════════════════════════════════════════════════════════════════════════
# Render orchestrator
# ═════════════════════════════════════════════════════════════════════════════

def render(run_name: str, results_dir: Path, samplesheet: Optional[Path],
           branding: Path, output: Path) -> str:
    rows = collect(results_dir, samplesheet)
    org  = load_org_info(branding)
    logo = encode_logo_b64(branding)
    epicandi_logo = file_to_b64(branding / "epicandi.png")
    batch_colors: dict[str, str] = {}

    quast_rows = collect_quast(results_dir)
    busco_rows = collect_busco(results_dir)

    # Secondary HTML for the full QC table.
    secondary = output.parent / f"{output.stem}_fastq_qc_full.html"
    render_secondary_qc_html(rows, secondary)

    secondary_asm = output.parent / f"{output.stem}_assembly_qc_full.html"
    render_secondary_assembly_html(quast_rows, busco_rows, secondary_asm)

    # Excel companion.
    excel = output.with_suffix(".xlsx")
    try:
        write_excel_companion(rows, quast_rows, busco_rows, results_dir, excel)
    except Exception as e:
        print(f"[generate_epicandi_report] Excel skipped: {e}", file=sys.stderr)

    sections = [
        render_banner(run_name, org, logo, epicandi_logo),
        render_summary_cards(rows),
        render_species_batch_section(rows, batch_colors),
        render_identification_chart(rows, batch_colors),
        render_fastq_qc(rows, batch_colors, secondary),
        render_assembly_qc(quast_rows, busco_rows, secondary_asm),
        render_amr_cnv(results_dir, rows, batch_colors),
        render_snp_section(results_dir, rows, batch_colors,
                           export_dir=output.parent),
    ]

    # Batch legend (if any batches with assigned colours)
    legend_html = ""
    if batch_colors:
        legend_html = ('<section class="section"><h3>Batch legend</h3>'
                       '<div class="legend">'
                       + "".join(f'<span><span class="sw" style="background:{c};"></span>{b}</span>'
                                 for b, c in batch_colors.items())
                       + '</div></section>')

    footer = ('<div class="footer">epicandi &middot; '
              f'{org.get("name","")} &middot; {datetime.now().strftime("%Y-%m-%d")}'
              f' &middot; <a href="{excel.name}">Excel companion</a></div>')

    # Plotly.js — preferentemente embebido inline desde branding/vendor para que
    # el reporte funcione offline (file:// o entornos sin red).  Si no está
    # disponible localmente, se recurre al CDN.
    # Embed Plotly + D3 inline so the report is a single self-contained file
    # (no sibling .js to forget when transferring via SFTP).  Both libraries
    # are minified — Plotly ~4.5 MB, D3 ~280 KB.
    plotly_cdn = ""
    for lib_name in ("plotly-2.35.2.min.js", "d3.v7.min.js"):
        src = branding / "vendor" / lib_name
        if src.exists():
            js_text = src.read_text(encoding="utf-8", errors="ignore")
            # Defensive: minified bundles never contain </script> but if any
            # third-party bundle ever does, splitting it avoids closing the
            # parent tag prematurely.
            js_text = js_text.replace("</script>", "<\\/script>")
            plotly_cdn += f'<script>{js_text}</script>'
        else:
            url = ("https://cdn.plot.ly/plotly-2.35.2.min.js"
                   if "plotly" in lib_name else "https://d3js.org/d3.v7.min.js")
            plotly_cdn += f'<script src="{url}" charset="utf-8"></script>'

    # After every chart has been created, force a resize for all Plotly divs.
    # This catches the rare case where a container's width was 0 at the moment
    # of newPlot (e.g. due to layout settling later), which results in invisible
    # bars / clipped axes.  Triggered on both load + a 500 ms safety net.
    plotly_cdn += (
        "<script>"
        "function _resizeAllPlots(){ "
        "  document.querySelectorAll('[id^=plotly_div_]').forEach(function(el){ "
        "    try { Plotly.Plots.resize(el); } catch(_){}; "
        "  });"
        "}"
        "window.addEventListener('load',  function(){ setTimeout(_resizeAllPlots, 100); });"
        "window.addEventListener('resize', _resizeAllPlots);"
        "setTimeout(_resizeAllPlots, 500);"
        "</script>"
    )

    # Visible diagnostic banner — shows library load status and surfaces any
    # uncaught JS error so the user can see problems without opening DevTools.
    plotly_cdn += (
        "<script>"
        "window.addEventListener('error', function(e){"
        "  var b = document.createElement('div');"
        "  b.style.cssText = 'position:sticky;top:0;z-index:9999;"
        "    background:#FDE2E2;color:#900;padding:8px 12px;"
        "    border-bottom:2px solid #c00;font:12px/1.4 monospace;';"
        "  b.textContent = '[JS error] ' + (e.message || e.error) + "
        "    '  ·  ' + (e.filename||'') + ':' + (e.lineno||'');"
        "  document.body && document.body.prepend(b);"
        "});"
        "window.addEventListener('load', function(){"
        "  var ok = (typeof Plotly !== 'undefined') && (typeof d3 !== 'undefined');"
        "  if (!ok) {"
        "    var b = document.createElement('div');"
        "    b.style.cssText = 'position:sticky;top:0;z-index:9999;"
        "      background:#FFF8E0;color:#a50;padding:8px 12px;"
        "      border-bottom:2px solid #d80;font:12px/1.4 monospace;';"
        "    b.textContent = '[libs] Plotly: ' + (typeof Plotly !== 'undefined' ? 'OK' : 'MISSING') + "
        "      '  ·  D3: ' + (typeof d3 !== 'undefined' ? 'OK' : 'MISSING');"
        "    document.body.prepend(b);"
        "  }"
        "});"
        "</script>"
    )

    return ("<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"utf-8\">"
            f"<title>epicandi report — {run_name}</title>"
            f"{plotly_cdn}"
            f"<style>{CSS}</style></head><body>"
            f"<div class=\"wrap\">{''.join(sections)}{legend_html}{footer}</div>"
            "</body></html>")


# ═════════════════════════════════════════════════════════════════════════════
# CLI
# ═════════════════════════════════════════════════════════════════════════════

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", required=True, type=Path)
    ap.add_argument("--run-name",    required=True)
    ap.add_argument("--samplesheet", type=Path, default=None)
    ap.add_argument("--branding",    type=Path, default=Path("branding"))
    ap.add_argument("--output",      required=True, type=Path)
    args = ap.parse_args()

    if not args.results_dir.exists():
        sys.exit(f"[generate_epicandi_report] results-dir not found: {args.results_dir}")

    html = render(args.run_name, args.results_dir, args.samplesheet,
                  args.branding, args.output)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(html, encoding="utf-8")
    print(f"[generate_epicandi_report] wrote {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

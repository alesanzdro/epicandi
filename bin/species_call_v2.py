#!/usr/bin/env python3
"""
species_call_v2.py — combine Sylph + READ_SCREEN (illumina + nanopore)
+ AuriClass + manifest lookup to emit a single species_call.tsv per sample.

Design (post-bugfix 2026-05-18):
  1. SYLPH is authoritative for species_assigned (most sensitive taxonomic caller).
     Aggregated AT SPECIES LEVEL (sum all slug Taxonomic_abundance per ncbi_name).
  2. READ_SCREEN reports per-read purity (the % reads metric the boss asked for).
     Also aggregated AT SPECIES LEVEL (sum reads_mapped per ncbi_name) — handles
     cross-clade homology (e.g. C. auris reads spread across cladeI/II/III slugs).
  3. AURICLASS only refines clade (does not change species_assigned).
  4. Classification combines both signals (see decision table below).

Output columns:
  sample_id, species_assigned (ncbi_name), species_clinical (legacy_name),
  classification {high_conf, mixed_culture, low_conf, novel_or_contamination},
  sylph_ani, sylph_tax_pct,
  readscreen_top_pct, readscreen_other_pct, readscreen_unmapped_pct,
  readscreen_total_reads,
  reference_slug (top Sylph hit slug — for cross-check vs samplesheet.reference_tag),
  ploidy, ploidy_conf, clade
"""
import argparse
import csv
import re
import sys
from pathlib import Path


def read_manifest(path):
    """Return dict keyed by any known slug alias (tag, sylph_slug, or legacy
    `slug` column) -> record. Multiple keys point to the same row so that
    SCREEN_TALLY output (uses `tag` as the per-ref header prefix) and Sylph
    output (uses `sylph_slug`) can both be resolved against the manifest.
    """
    out = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            for col in ("tag", "sylph_slug", "slug"):
                key = (row.get(col) or "").strip()
                if key:
                    out[key] = row
    return out


def slug_from_path(path_str):
    """Extract slug from any path string — looks for the LAST path component
    that matches a known slug pattern. Strategy: try every basename and parent
    dir name; the caller has the manifest to validate.

    Sylph DB might emit <sylph_db_dir>/cladeI_B8441.fna.gz
    or <refs_dir>/auris/cladeI_B8441/genome.fna.gz.
    Both have 'cladeI_B8441' as a recognizable component."""
    parts = Path(path_str).parts
    candidates = []
    for p in parts:
        # Strip extensions
        base = p.replace('.fna.gz', '').replace('.fasta.gz', '').replace('.fasta', '').replace('.fna', '')
        candidates.append(base)
    return candidates


# Slug aliases — Sylph DB was sometimes built with a different strain reference
# than the panel/manifest. These map Sylph-emitted slugs to manifest slugs
# (same species, different type-strain).
SLUG_ALIASES = {
    "Ckrusei_CBS573": "Ckrusei_ATCC_6258",   # Pichia kudriavzevii — same species
}


def find_slug_in_manifest(path_str, manifest_slugs):
    """Return the slug if any path component matches a known manifest slug.
    Applies SLUG_ALIASES to handle DB/panel strain mismatches."""
    for cand in slug_from_path(path_str):
        # Direct match
        if cand in manifest_slugs:
            return cand
        # Alias match
        if cand in SLUG_ALIASES and SLUG_ALIASES[cand] in manifest_slugs:
            return SLUG_ALIASES[cand]
    return None


def parse_sylph_profile(path):
    """Sylph profile output is TSV."""
    if not Path(path).exists():
        return []
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip().split("\t")
        for line in fh:
            parts = line.rstrip().split("\t")
            d = dict(zip(header, parts))
            rows.append(d)
    return rows


def parse_readscreen(path):
    """READ_SCREEN tally TSV: sample, reference, reads_mapped.
    Reference is prefixed with 'EpiCandi_' — strip it to get the slug."""
    if not path or not Path(path).exists():
        return {}, 0
    by_slug = {}
    total = 0
    with open(path) as fh:
        header = fh.readline().rstrip().split("\t")
        for line in fh:
            parts = line.rstrip().split("\t")
            d = dict(zip(header, parts))
            ref = d.get("reference", "")
            cnt = int(d.get("reads_mapped", 0) or 0)
            slug = re.sub(r'^EpiCandi_', '', ref)
            if slug.lower() == "unmapped":
                by_slug["__unmapped__"] = cnt
            else:
                by_slug[slug] = by_slug.get(slug, 0) + cnt
            total += cnt
    return by_slug, total


def parse_auriclass(path):
    if not path or not Path(path).exists():
        return None
    with open(path) as fh:
        header = fh.readline().rstrip().split("\t")
        line = fh.readline().rstrip()
        if not line:
            return None
        parts = line.split("\t")
        d = dict(zip(header, parts))
    raw = d.get("Clade") or d.get("clade") or (parts[1] if len(parts) > 1 else "")
    clade = raw.replace("Clade ", "").strip()
    return clade or None


def aggregate_by_species(by_slug, manifest):
    """Map slug counts/abundances to species-level (ncbi_name)."""
    by_species = {}
    species_top_slug = {}   # which slug had highest value per species (for reference_slug)
    species_top_val = {}
    for slug, val in by_slug.items():
        if slug == "__unmapped__":
            continue
        row = manifest.get(slug)
        if row is None:
            continue
        ncbi = row.get("ncbi_name", "Unknown")
        by_species[ncbi] = by_species.get(ncbi, 0) + val
        if val > species_top_val.get(ncbi, -1):
            species_top_val[ncbi] = val
            species_top_slug[ncbi] = slug
    return by_species, species_top_slug


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample_id", required=True)
    ap.add_argument("--sylph", required=True)
    ap.add_argument("--readscreen_illumina")
    ap.add_argument("--readscreen_nanopore")
    ap.add_argument("--auriclass")
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--ani_high", type=float, default=97.0)
    ap.add_argument("--ani_low", type=float, default=95.0)
    ap.add_argument("--readscreen_pct_top", type=float, default=85.0,
                    help="Min %% reads to top species for high_conf")
    ap.add_argument("--readscreen_pct_unmapped_max", type=float, default=50.0,
                    help="Max %% unmapped reads before flagging novel_or_contamination")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    manifest = read_manifest(args.manifest)
    slugs = set(manifest.keys())
    sylph_rows = parse_sylph_profile(args.sylph)

    # === SYLPH AGGREGATION ===
    sylph_by_slug = {}
    sylph_ani_by_slug = {}
    for r in sylph_rows:
        slug = find_slug_in_manifest(r.get("Genome_file", ""), slugs)
        if slug is None:
            continue
        tax_pct = float(r.get("Taxonomic_abundance", 0) or 0)
        ani = float(r.get("Adjusted_ANI", r.get("ANI", 0)) or 0)
        sylph_by_slug[slug] = sylph_by_slug.get(slug, 0) + tax_pct
        if ani > sylph_ani_by_slug.get(slug, -1):
            sylph_ani_by_slug[slug] = ani

    sylph_by_species, sylph_top_slug_per_species = aggregate_by_species(sylph_by_slug, manifest)

    if not sylph_by_species:
        # FALLBACK: try READ_SCREEN if Sylph emitted nothing. Sample may have
        # divergent strain that Sylph's k-mer DB couldn't ANI-confidently call,
        # but bowtie2 alignment to the panel still captured the species.
        rs_by_slug_fb = {}
        rs_total_fb = 0
        for p in (args.readscreen_illumina, args.readscreen_nanopore):
            if not p:
                continue
            bs, t = parse_readscreen(p)
            for slug, cnt in bs.items():
                rs_by_slug_fb[slug] = rs_by_slug_fb.get(slug, 0) + cnt
            rs_total_fb += t
        rs_unmapped_fb = rs_by_slug_fb.pop("__unmapped__", 0)
        rs_by_species_fb, top_slugs_fb = aggregate_by_species(rs_by_slug_fb, manifest)

        if rs_by_species_fb:
            top_sp = max(rs_by_species_fb, key=rs_by_species_fb.get)
            top_slug = top_slugs_fb[top_sp]
            top_row = manifest[top_slug]
            top_pct = 100.0 * rs_by_species_fb[top_sp] / rs_total_fb if rs_total_fb else 0
            unmap_pct = 100.0 * rs_unmapped_fb / rs_total_fb if rs_total_fb else 0
            # READ_SCREEN >= 60% to top → low_conf (Sylph could not confirm)
            classif = "low_conf" if top_pct >= 60 else "novel_or_contamination"
            print(f"[species_call_v2] {args.sample_id}: Sylph empty — "
                  f"READ_SCREEN fallback → {top_sp} ({classif}, "
                  f"rs_top_pct={top_pct:.1f}%, rs_unmapped={unmap_pct:.1f}%)",
                  file=sys.stderr)
            write_row_full(args, top_row, classif,
                           sylph_ani=0.0, sylph_tax_pct=0.0,
                           rs_top_pct=round(top_pct, 2),
                           rs_other_pct=round(100.0 - top_pct - unmap_pct, 2),
                           rs_unmapped_pct=round(unmap_pct, 2),
                           rs_total=rs_total_fb,
                           reference_slug=top_slug,
                           clade="")
            return

        write_unknown(args, reason="no Sylph hit AND no READ_SCREEN signal")
        return

    sylph_top_species = max(sylph_by_species, key=sylph_by_species.get)
    sylph_tax_pct = sylph_by_species[sylph_top_species]
    sylph_top_slug = sylph_top_slug_per_species[sylph_top_species]
    sylph_top_ani = sylph_ani_by_slug.get(sylph_top_slug, 0)

    # Manifest row for the top hit (drives reference_slug, ploidy, ploidy_conf, is_auris)
    top_row = manifest[sylph_top_slug]
    species_assigned = top_row.get("ncbi_name", "Unknown")
    species_clinical = top_row.get("legacy_name", "Unknown")
    reference_slug = sylph_top_slug
    ploidy = top_row.get("ploidy", "1")
    ploidy_conf = top_row.get("ploidy_conf", "baja")
    is_auris = top_row.get("is_auris", "false").lower() == "true"

    # === READ_SCREEN AGGREGATION ===
    rs_by_slug = {}
    rs_total = 0
    for p in (args.readscreen_illumina, args.readscreen_nanopore):
        if not p:
            continue
        bs, t = parse_readscreen(p)
        for slug, cnt in bs.items():
            rs_by_slug[slug] = rs_by_slug.get(slug, 0) + cnt
        rs_total += t

    rs_unmapped = rs_by_slug.pop("__unmapped__", 0) if "__unmapped__" in rs_by_slug else 0
    rs_by_species, _ = aggregate_by_species(rs_by_slug, manifest)
    rs_mapped_total = sum(rs_by_species.values())

    if rs_total > 0:
        rs_top_pct = round(100.0 * rs_by_species.get(species_assigned, 0) / rs_total, 2)
        rs_other_pct = round(100.0 * (rs_mapped_total - rs_by_species.get(species_assigned, 0)) / rs_total, 2)
        rs_unmapped_pct = round(100.0 * rs_unmapped / rs_total, 2)
    else:
        rs_top_pct = rs_other_pct = rs_unmapped_pct = 0.0

    # === CLASSIFICATION ===
    if rs_unmapped_pct > args.readscreen_pct_unmapped_max:
        classification = "novel_or_contamination"
    elif sylph_top_ani < args.ani_low:
        classification = "low_conf"
    elif sylph_top_ani >= args.ani_high and rs_top_pct >= args.readscreen_pct_top:
        classification = "high_conf"
    elif sylph_top_ani >= args.ani_low and (rs_other_pct >= 5 or rs_top_pct < 70):
        classification = "mixed_culture"
    elif sylph_top_ani >= args.ani_high:
        # Sylph confident but READ_SCREEN slightly under threshold
        # (common for closely-related clades like auris) → still high_conf
        classification = "high_conf"
    else:
        classification = "low_conf"

    # === CLADE (auris only) ===
    clade = ""
    if is_auris:
        ac = parse_auriclass(args.auriclass) if args.auriclass else None
        if ac:
            clade = ac
        else:
            slug = reference_slug
            if slug.startswith("clade"):
                clade = slug.split("_")[0].replace("clade", "")

    # === WRITE ===
    header = [
        "sample_id", "species_assigned", "species_clinical", "classification",
        "sylph_ani", "sylph_tax_pct",
        "readscreen_top_pct", "readscreen_other_pct", "readscreen_unmapped_pct",
        "readscreen_total_reads",
        "reference_slug", "ploidy", "ploidy_conf", "clade",
    ]
    row = [
        args.sample_id, species_assigned, species_clinical, classification,
        str(round(sylph_top_ani, 2)), str(round(sylph_tax_pct, 4)),
        str(rs_top_pct), str(rs_other_pct), str(rs_unmapped_pct),
        str(rs_total),
        reference_slug, ploidy, ploidy_conf, clade,
    ]
    with open(args.output, "w") as fh:
        fh.write("\t".join(header) + "\n")
        fh.write("\t".join(row) + "\n")

    print(f"[species_call_v2] {args.sample_id}: {species_assigned} "
          f"({classification}, ANI={sylph_top_ani}, sylph_tax={sylph_tax_pct:.2f}%, "
          f"rs_top={rs_top_pct}%, rs_other={rs_other_pct}%, rs_unmap={rs_unmapped_pct}%, "
          f"ploidy={ploidy}{', clade '+clade if clade else ''})",
          file=sys.stderr)


def write_row_full(args, manifest_row, classification, sylph_ani, sylph_tax_pct,
                   rs_top_pct, rs_other_pct, rs_unmapped_pct, rs_total,
                   reference_slug, clade):
    species_assigned = manifest_row.get("ncbi_name", "Unknown")
    species_clinical = manifest_row.get("legacy_name", "Unknown")
    ploidy           = manifest_row.get("ploidy", "1")
    ploidy_conf      = manifest_row.get("ploidy_conf", "baja")
    header = [
        "sample_id", "species_assigned", "species_clinical", "classification",
        "sylph_ani", "sylph_tax_pct",
        "readscreen_top_pct", "readscreen_other_pct", "readscreen_unmapped_pct",
        "readscreen_total_reads",
        "reference_slug", "ploidy", "ploidy_conf", "clade",
    ]
    row = [
        args.sample_id, species_assigned, species_clinical, classification,
        str(round(sylph_ani, 2)), str(round(sylph_tax_pct, 4)),
        str(rs_top_pct), str(rs_other_pct), str(rs_unmapped_pct),
        str(rs_total),
        reference_slug, ploidy, ploidy_conf, clade,
    ]
    with open(args.output, "w") as fh:
        fh.write("\t".join(header) + "\n")
        fh.write("\t".join(row) + "\n")


def write_unknown(args, reason="no signal"):
    header = [
        "sample_id", "species_assigned", "species_clinical", "classification",
        "sylph_ani", "sylph_tax_pct",
        "readscreen_top_pct", "readscreen_other_pct", "readscreen_unmapped_pct",
        "readscreen_total_reads",
        "reference_slug", "ploidy", "ploidy_conf", "clade",
    ]
    row = [
        args.sample_id, "Unknown", "Unknown", "novel_or_contamination",
        "0", "0",
        "0", "0", "0",
        "0",
        "", "1", "baja", "",
    ]
    with open(args.output, "w") as fh:
        fh.write("\t".join(header) + "\n")
        fh.write("\t".join(row) + "\n")
    print(f"[species_call_v2] {args.sample_id}: Unknown ({reason})", file=sys.stderr)


if __name__ == "__main__":
    main()

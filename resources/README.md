# `resources/` — external data assets

This directory holds large external assets (references, databases, models)
that the pipeline consumes but does not version-control in git. The
`.gitignore` excludes everything except this README.

## Layout

```
resources/
├── references/                              # one dir per reference_tag
│   └── <tag>/
│       ├── reference.fasta                  # genome FASTA
│       ├── reference.fasta.fai              # samtools index (auto-built if absent)
│       ├── reference.gff3                   # funannotate GFF (for coord lookups)
│       └── mask.bed                         # rRNA + repeats + low-complexity
│                                            # → consumed as --exclude-intervals by HC
├── databases/
│   ├── fungamr/
│   │   └── FungAMR_070425_epicandi.csv      # FungAMR cross-reference for build_mutation_catalog.py
│   ├── sylph/
│   │   └── sketches/
│   │       ├── epicandi.syldb               # sylph reference DB sketch
│   │       └── ...
│   ├── clair3_models/
│   │   └── r1041_e82_400bps_sup_v500/       # PyTorch model (.pt files)
│   └── busco_downloads/
│       └── lineage/
│           └── saccharomycetes_odb12/
```

## Bootstrap on a new host

Re-create the symlinks/copies as below — paths shown match the FISABIO Heimdal
layout used by Alex; on a different host adjust accordingly.

```bash
mkdir -p resources/{references,databases/{fungamr,sylph,clair3_models,busco_downloads}}

# References (13 strains) — symlink fasta + .fai + gff3 + mask.bed
for tag in Calbicans_SC5314 CaurisI_B8441 CaurisII_B11220 CaurisIII_B11221 \
           CaurisIV_B11245 CaurisV_IFRC2087 CaurisVI_F3485 \
           Cdubliniensis_CD36 Clusitaniae_P1 Cparapsilosis_CDC317 \
           Ctropicalis_MYA3404 Nglabratus_CBS138 Pkudriavzevii_CBS573T; do
    mkdir -p resources/references/${tag}
    ln -sf /home/asanzc/0epicandi/v1.7/00_INPUT/references/${tag}/reference.fasta     resources/references/${tag}/
    ln -sf /home/asanzc/0epicandi/v1.7/00_INPUT/references/${tag}/reference.fasta.fai resources/references/${tag}/ 2>/dev/null
    ln -sf /home/asanzc/0epicandi/v1.7/02_OUTPUT/${tag}/reference.gff3                resources/references/${tag}/
    ln -sf /home/asanzc/0epicandi/v1.7/02_OUTPUT/${tag}/mask.bed                       resources/references/${tag}/
done

# DBs (symlinks; copy if portability matters more than disk)
ln -sfn /home/asanzc/epicandi-nf/resources/FungAMR_070425_epicandi.csv  resources/databases/fungamr/
ln -sfn /home/asanzc/epicandi_databases/sylph/sketches                  resources/databases/sylph/sketches
ln -sfn /home/asanzc/epicandi/resources/clair3_models/r1041_e82_400bps_sup_v500  resources/databases/clair3_models/
ln -sfn /home/asanzc/epicandi/resources/busco_downloads                 resources/databases/busco_downloads/lineage
```

## Why symlinks (not copies)?

The references are ~1.5 GB total (13 × ~12 MB fasta + GFF + masks); the sylph
DB is ~50 MB; busco saccharomycetes_odb12 is ~500 MB; clair3 model ~20 MB.
Total ≈ 750 MB. Symlinks keep the repo small and let the same data serve
multiple branches/clones on this host. For a truly portable bundle, replace
the `ln -sf` with `cp` and ship the directory.

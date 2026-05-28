# Test data

## Smoke test (Illumina, single sample)

The `-profile test` smoke test is wired in `conf/test.config` to read the
samplesheet at `assets/samplesheet.csv`, which points to:

```
/almacenamiento/asanzc/SRA_EXP_CANDIDA/03_fastq/Cauris_B11220/SRX26245890/
  SRR30847274_1.fastq.gz   (~480 MB)
  SRR30847274_2.fastq.gz   (~475 MB)
```

This is a *C. auris* B11220 (clade II) Illumina paired-end sample.  It lives
on `/almacenamiento/` (the institutional NAS) and is accessible on any host
that mounts that volume.

If you run the pipeline from a different host:

1. Download SRR30847274 from the NCBI SRA (~950 MB), or
2. Edit `assets/samplesheet.csv` to point to your local copy.

### Smoke test runtime

≈ 15-25 min on Heimdal/Gru (16 cores) with the conda envs already cached.

The smoke is single-sample, so the cohort-level steps automatically skip
with a warning:

- Joint genotyping (CombineGVCFs / GenotypeGVCFs) needs `params.cohort_min_samples` (default 3) samples in the same cohort.
- Phylogeny + UPGMA network need ≥ 3 samples.

The smoke does cover: QC, identification, assembly, post-asm QC, AMR,
per-sample variant calling, CNV, and the HTML report.

## Nanopore subsampled test (TBD)

A Nanopore-only mini sample (~30× *C. auris*, ~370 MB FASTQ) is on the
roadmap.  Suggested source: subsample an existing run from
`/almacenamiento/asanzc/SRA_EXP_CANDIDA/` using
`filtlong --target_bases <gsize*30>` and place it under this directory.

Once added, create a `samplesheet_nanopore_test.csv` and a
`conf/test_nanopore.config` profile pointing to it.

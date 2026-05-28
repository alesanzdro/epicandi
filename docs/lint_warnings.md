# Lint warnings — accepted

`nf-core pipelines lint` for `epicandi/epicandi` reports five warnings that are intentionally not resolved. Documented here so the next maintainer doesn't waste cycles re-investigating.

## 1. `pipeline_name_conventions` — dash in pipeline name

```
pipeline_name_conventions: Naming does not adhere to nf-core conventions:
Contains non alphanumeric characters
```

Pipeline name is `epicandi-nf3` (host directory `/home/asanzc/epicandi-nf3/`). nf-core convention is lowercase alphanumeric without separator. The dash was chosen at CHECKPOINT A to disambiguate from prior `epicandi-nf` and `epicandi-nf2` attempts on the same HPC. Renaming would break the manifest/homepage URL and the directory layout the user explicitly picked.

Mitigation: workflow include aliases inside `main.nf` use underscores (`EPICANDI`, `EPICANDI_EPICANDI`) so the Nextflow parser never sees the dash.

## 2. `qc_gate main_nf_include` — subworkflow includes < 2 modules

```
Subworkflow includes less than two modules
```

`subworkflows/local/qc_gate/` only wraps the `QC_GATE` module, but it adds non-trivial logic (per-sample join of clean illumina + nanopore reads with `NO_FILE_*` placeholders, parse of the PASS/WARN/FAIL flag, branching of the samplesheet to drop FAIL samples). Inlining that join into the top-level workflow would clutter `workflows/epicandi.nf`. nf-core's "≥2 modules per subworkflow" heuristic is a stylistic suggestion, not a correctness requirement.

## 3. `qc_gate meta_name` — workflow alias suffix

```
Conflicting workflow name between meta.yml (`qc_gate`) and main.nf (`QC_GATE_SWF`)
```

The workflow is named `QC_GATE_SWF` because the module it wraps is also called `QC_GATE` — Nextflow allows same-named entities in different namespaces but side-by-side reading is confusing. Renaming the dir to `qc_gate_swf/` would lose the alignment with the module's `qc/gate/` path; renaming the workflow to `QC_GATE` would force `include {QC_GATE as QC_GATE_SUB}` aliases at every call site.

## 4. `subsample meta_name` — workflow alias suffix

```
Conflicting workflow name between meta.yml (`subsample`) and main.nf (`SUBSAMPLE_SWF`)
```

Same pattern as `qc_gate`: the local `SUBSAMPLE` *module* (`modules/local/subsample/main.nf`) is wrapped by a `SUBSAMPLE_SWF` *subworkflow* (`subworkflows/local/subsample/main.nf`). Accepted for the same reason.

## 5. `utils_nfcore_epicandi_pipeline meta_name` — template helper

```
Conflicting workflow name between meta.yml (`utils_nfcore_epicandi_pipeline`)
and main.nf (`PIPELINE_INITIALISATION`)
```

`subworkflows/local/utils_nfcore_<pipeline>_pipeline/` is a template artefact shipped with `nf-core pipelines create`. It bundles two workflows in the same `main.nf` (`PIPELINE_INITIALISATION` + `PIPELINE_COMPLETION`) plus helper functions, which inherently violates the "one workflow per meta.yml, name = workflow" rule.

## Auto-skipped tests (in `.nf-core.yml`)

| Test | Reason |
|---|---|
| `container_configs` | `nf-core 4.0.2` lint cannot parse `nextflow inspect` output that Nextflow 25.x prefixes with `[PIPELINE]` / `[WORKDIR]` banner lines. Re-enable once upstream lint strips those banners. |
| `subworkflow_changes: subworkflows/local/utils_nfcore_epicandi_pipeline` | Template helper, pipeline-specific by design. |
| `merge_markers: modules/nf-core/bowtie2/align/tests/main.nf.test.snap` | Upstream test snapshot embeds 7+ consecutive `<` in a SAM quality string; the merge-marker regex misfires. |
| `files_unchanged: [logos, CODE_OF_CONDUCT.md, …]` | We do not maintain nf-core logos / CoC / awstest workflows. |
| `files_exist: [logos, awstest, …, conf/igenomes*.config]` | Not applicable to a non-nf-core, non-iGenomes pipeline. |
| `nextflow_config: manifest.name`, `manifest.homePage` | `epicandi/epicandi` is not a `nf-core/*` pipeline; lint check assumes the nf-core namespace. |

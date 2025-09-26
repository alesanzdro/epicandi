# Changelog

All notable changes to this project will be documented in this file.

The format is based on Keep a Changelog (https://keepachangelog.com/en/1.1.0/) 
and this project aims to adhere to Semantic Versioning (https://semver.org/).

## [Unreleased]
- Placeholder for upcoming changes.

## - 2025-09-26
### Added
- Configurable phylogenetic pipeline with selectable methodologies (none, triage, orthologs, all).
- Consolidated phylogenetic alignment processing script integrating selection, fixing and filtering steps.
- Auxiliary scripts for ANI distance conversion and Mash matrix transformation.
- Core ortholog selection and alignment quality processing helpers.
- Documentation covering IQ-TREE issues and full phylogenetic system update.

### Changed
- Phylogenetic activation mechanism moved to configuration key (`phylogenetics.methodology`).

### Removed
- Legacy secondary Snakefile variant in favor of unified configurable approach.
- Residual AI reference files and cached bytecode noise.

### Breaking Changes
- Phylogenies are now controlled exclusively through the configuration file instead of implicit rule invocation.

## - 2025-09-26
### Added
- Comprehensive dual phylogenetic analysis system (ortholog-based + distance-based) prior to consolidation.

### Changed
- Iterative refinement of pipeline structure preparing for configuration-driven control.

## - 2025-09-23
### Changed
- Major modernization and standardization refactor across rules and scripts.
- Modularization of assembly workflow with enhanced strategies.

## - 2025-09-22
### Added
- Modular assembly pipeline structure introducing clearer stage separation.

## - 2025-09-20
### Added
- Resistance analysis script for characterization layer.
- Initial full pipeline commit including primary Snakefile, configuration, and environment definitions.

### Documentation
- README updates and initial contributor guidance.

### Chore
- Initial .gitignore to manage environment and bytecode artifacts.

## [Deprecated / Historical]
- Temporary documentation file for AI assistance (removed in later cleanup).

---

### Version Derivation Notes
Version numbers are retroactively assigned based on commit chronology and scope magnitude:
- 0.1.x: Initial functional baseline and incremental scripts.
- 0.2.0: Structural assembly modularization.
- 0.3.0: Broad refactor and standardization.
- 0.4.0: Dual phylogeny system introduction.
- 0.5.0: Configuration-driven phylogeny with breaking change.

### Upgrade Guidance
- Upgrading to 0.5.0 requires setting `phylogenetics.methodology` in `config.yaml`.
- Prior versions relied on implicit rule availability; this is no longer valid.

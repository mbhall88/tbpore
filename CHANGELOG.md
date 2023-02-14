# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]

## [0.5.0]

### Added

- Read stats (subsampled reads) from [`nanoq`](https://github.com/esteinig/nanoq) in `tbpore process`. The stats file will be in the output directory.

## [0.4.1]

### Changed

- Removed 'subsampled' from the VCF path [[#49][49]]

## [0.4.0]

### Added

- `download` subcommand to download and validate the decontamination database
- sorting of the decontaminated fastq file to ensure subsampling is reproducible regardless of whether a combined fastq or directory of fastqs is provided [[#48][48]]

### Changed

- Default (expected) location of the decontamination database is now `${HOME}/.tbpore/decontamination_db/remove_contam.map-ont.mmi`
- Variant filtering params
  * minimum depth changed from 0 to 5
  * minimum variant quality (QUAL) score from 85 to 25
  * minimum mapping quality from 0 to 30

### Removed

- Usage of deprecated `iteritems` function from pandas. This will remove an annoying deprecated warning when running `tbpore cluster`

## [0.3.2]

### Changed

- When `--name` is not given, take name to be the filename minus the fastq (and optional gz) suffix. Previously, we took everything before the first `.` [[#45][45]]

## [0.3.1]

### Fixed

- `data/` and `external_scripts/` directories were not having their contents included in site-packages when installing from sdist

## [0.3.0]

### Added

- Ability to specify the path to the cache directory (`--cache`) [[#43][43]]

### Changed

- Default cache dir is now `${HOME}/.cache` [[#43][43]]
- Default output directory for all subcommands is now the current directory

## [0.2.0]

### Added

- Ability to specify a different path for the decontamination database (`--db`) and metadata file (`--metadata`) [[#34]][34]

### Changed

- Update mykrobe to v0.12 (WHO catalogue)

## [0.1.1]

### Fixed

- Added missing, required files and directories, `.config.yaml`, `data/`, `external_scripts/`, and `CHANGELOG.md` to distribution

## [0.1.0]

### Added

- First release - so everything you see is new!

[unreleased]: https://github.com/mbhall88/tbpore/compare/0.5.0...HEAD
[0.5.0]: https://github.com/mbhall88/tbpore/compare/0.4.1...0.5.0
[0.4.1]: https://github.com/mbhall88/tbpore/compare/0.4.0...0.4.1
[0.4.0]: https://github.com/mbhall88/tbpore/compare/0.3.2...0.4.0
[0.3.2]: https://github.com/mbhall88/tbpore/compare/0.3.1...0.3.2
[0.3.1]: https://github.com/mbhall88/tbpore/compare/0.3.0...0.3.1
[0.3.0]: https://github.com/mbhall88/tbpore/compare/0.2.0...0.3.0
[0.2.0]: https://github.com/mbhall88/tbpore/compare/0.1.1...0.2.0
[0.1.1]: https://github.com/mbhall88/tbpore/compare/0.1.0...0.1.1
[0.1.0]: https://github.com/mbhall88/tbpore/releases/tag/0.1.0
[34]: https://github.com/mbhall88/tbpore/issues/34
[43]: https://github.com/mbhall88/tbpore/issues/43
[45]: https://github.com/mbhall88/tbpore/issues/45
[48]: https://github.com/mbhall88/tbpore/issues/48
[49]: https://github.com/mbhall88/tbpore/issues/49

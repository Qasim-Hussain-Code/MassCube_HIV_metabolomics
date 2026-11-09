# data/

Nothing in this directory is tracked by git except this file. Every file
here is fetched again by the pipeline and checked against a checksum that
the data source publishes, so a fresh clone recreates it exactly.

| path | written by | what it holds | checked against |
|---|---|---|---|
| `metadata/` | `scripts/02_download_metadata.sh` | the MTBLS9264 ISA-Tab (`i_`, `s_` and two `a_` files), the authors' two metabolite assignment files (`m_..._maf.tsv`), the published checksum lists and the MetaboLights study record | sha256 in `HASHES/metadata_sha256.json` |
| `raw/pos/`, `raw/neg/` | `scripts/04_process_mode.sh` | one Thermo `.raw` file at a time, deleted as soon as it has been converted and processed | sha256 in `HASHES/data_sha256.json` |
| `db/` | the annotation stage | the MassCube MS/MS database pickle from Zenodo, and the indexed subset of the GNPS Drug Library | md5 in the Zenodo record, recorded in `config/tool_versions.tsv` |

The raw files are about 300 MB each and 72,715,682,516 bytes in all
(`results/mtbls9264_study.tsv`), more than the 10 GB this project allows
itself, which is why they never accumulate here: each one is downloaded,
verified, converted to mzML, processed by MassCube and deleted before the
next is fetched. The converted mzML files live briefly under `work/` and
are deleted the same way. What stays is listed per file in
`results/conversion_manifest.tsv`: the raw and mzML sha256, sizes, scan
counts, polarity, instrument and acquisition time, which are enough to
reproduce and check each conversion.

The authors' metabolite assignment files are downloaded with the rest of
the metadata but are not read until the primary analysis is finished
(`preregistration.md`, section 1).

The data are distributed by MetaboLights under the EMBL-EBI Terms of Use
(the study's licence field, revision 2 of 8 September 2025).

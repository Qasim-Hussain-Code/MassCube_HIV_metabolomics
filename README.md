# MassCube_HIV_metabolomics

An analysis of MetaboLights MTBLS9264, untargeted LC-MS of stool from people living with HIV and healthy controls, processed from the raw spectra with MassCube. The pipeline and the results will be added here as each stage is finished.

## Status

As of 12 November 2026 the pipeline is built and pre-registered, and two of the 229 raw files have been processed. Every stage from configuration to figures is written and committed. Stages 00 to 06 have run on real data: the environment, the study metadata, the sample tables, and the per-file path (download, checksum, conversion, mzML checks, MassCube) for the positive-mode pooled QCs QC1 and QC11, with the pilot rules and the alignment tested on those two files. Stages 07 to 12 (annotation, filters, statistics, the published claims, the exploratory analyses and the figures) have been run only on simulated tables in the same format, to test the code; no result in this repository yet compares people living with HIV (PLWH) with healthy controls (HC).

The rest of the run is waiting on the build machine's disk. On 12 November, after Windows had grown its pagefile, drive C: no longer had room for one mode's per-file MassCube outputs above the 2 GB of free space the pre-registration keeps, and the EBI download rate had fallen back from megabytes to a few kilobytes per second. `bash run_all.sh` resumes where it stopped once space allows, and every stage refuses to start a download that would break the floor or the 10 GB ceiling.

## Pre-registration

[`preregistration.md`](preregistration.md) was written after reading the paper, its supplement and the ISA-Tab metadata, and committed on 9 November 2026 (`6d79eb5`, tag `preregistration`) before any spectrum was downloaded. It fixes the question, every processing parameter and QC rule, the primary test (limma on log2 peak heights, Benjamini-Hochberg within each ion mode, q < 0.05), the sensitivity analyses, the annotation levels, the antiretroviral positive control, the verdict rules for each published claim, and the storage plan. It is not edited; departures are listed under [Deviations from the pre-registration](#deviations-from-the-pre-registration). The authors' two metabolite assignment files were downloaded with the metadata and are opened only by the replication stage, after the primary analysis.

## Data

| | MTBLS9264 |
|---|---|
| Organism, material | Homo sapiens, stool |
| Publication | Tian X, Xie Y, Yu L, Yao P, Dong M, Jin C, Wu N. Microbiology Spectrum 2024, 12(11):e00238-24 |
| Accession | MetaboLights MTBLS9264, revision 2 of 8 September 2025 |
| Design | 70 PLWH (32 immune non-responders, 38 immune responders) and 34 HC, cross-sectional, Zhejiang, November 2020 to December 2022 |
| Instrument (from the raw files) | Thermo Q Exactive HF-X Orbitrap, Vanquish UHPLC, positive and negative ESI as separate runs |
| Raw files | 229, 72,715,682,516 bytes: 115 positive (38,138,325,295 bytes in the assay table plus `QC11.raw`, 350,562,367 bytes) and 114 negative (34,226,794,854 bytes) |
| Pooled QCs | 11 positive, 10 negative; no blanks deposited |
| Per-sample covariates | none deposited (no age, sex, ART regimen, CD4 count or viral load) |

The counts come from `results/mtbls9264_study.tsv` and `results/sample_counts.tsv`. Every raw file has a sha256 in the deposit's `HASHES/data_sha256.json`, and each one is checked against it before conversion.

Several problems in the deposit were found before any spectrum was read, and the first two decide how the sample tables are built.

The positive assay table's Sample Name column is out of step with every other identifier column in all 114 rows. Extract Name, MS Assay Name, Derived Spectral Data File and Normalization Name all equal the raw file's name, while Sample Name holds the same 114 names in another order, the ten QCs first; read literally it would put 40 files in the wrong group (`results/isa_assay_consistency.tsv`). The negative table agrees in every column. Files are linked to participants through Extract Name, which also agrees with the group prefix in every file name, and the antiretroviral control is computed under both links as a check.

`QC11.raw` is on the FTP with a published sha256 but in neither assay table, and the study table lists QC1 to QC10 only. Its scans are all of positive polarity (`results/conversion_manifest.tsv`), so it is kept as the eleventh positive-mode QC, as the pre-registration allows, and a sensitivity analysis drops it.

The method descriptions disagree with each other and with the files. The paper gives a Q Exactive HF-X (which it calls a column), 60,000 resolution, 8 scans per cycle, a 150 mm HSS T3 column at 0.25 mL/min, and gradients of 14 and 17 min, one step of which reads "1,420 min". The ISA protocol gives a Q Exactive, 70,000 resolution, 10 scans per cycle, a 100 mm column at 0.3 mL/min and a 12 min gradient starting at 8 per cent organic, and its data transformation text describes GC-MS data processed in MS-DIAL against a GC library, which cannot be the processing of these files. The two positive QCs read so far were acquired on a "Q Exactive HF-X Orbitrap" with the instrument method `20201222-16.5min-pos-vanquish(A1-0.1%FA+B1-0.1%FA-ACN)`, and their last scans are at 16.50 min (`results/conversion_manifest.tsv`): the instrument is the paper's, and the run is 16.5 min, which matches neither text.

Smaller points. The MetaboLights record reports a study size of 87,675,355,136 bytes, while the files on the FTP add up to 72,715,682,516 (`results/mtbls9264_study.tsv`). Table S1 gives the PLWH age as 42.48 years with a standard deviation of 2.17 against 11.59 for the controls; for 70 adults that is more likely a standard error, but the supplement does not say. The ISA extraction protocol adds 2-amino-3-(2-chlorophenyl)propionic acid at 4 ppm to every sample, which this pipeline uses as a process internal standard. The study's licence field is the EMBL-EBI Terms of Use.

## Pipeline

Numbered stages, run in order by `run_all.sh`. Each sources `project.conf`, skips itself when its stamp in `logs/` exists and says so, and writes its elapsed time, peak resident memory and the project's size on disk before and after to `logs/<stage>.resources.tsv`.

**Configure and install.** `00_configure.sh` detects or accepts threads, RAM and disk and writes `project.conf`. `01_install.sh` builds one conda environment inside the project (`env/`), from conda-forge only because bioconda has no Windows builds: Python 3.11 with MassCube 1.2.22 from PyPI (which pins ms_entropy 1.3.4, the version the MassCube database was built with), R 4.6.2 with limma, ropls, vegan and data.table installed by BiocManager inside that R, ThermoRawFileParser 1.4.5 from its GitHub release (its zip's sha256 is recorded at first retrieval in `config/external_checksums.tsv`), and ShellCheck. Every cache is kept inside the project and emptied afterwards, so that the disk budget measures it. The versions in use are in `config/tool_versions.tsv`, `config/python_packages.tsv`, `config/r_packages.tsv` and `config/environment.lock.yml`. Two things did not work as first written: MSYS `du` took more than half an hour over the environment's files with the virus scanner running, so sizes are read from directory listings by `scripts/lib/disk_usage.py`; and conda's cleaning keeps the extracted packages an environment uses, which on one drive are hard links into `env/` and were being counted twice, so the install deletes them.

**Metadata and sample tables.** `02_download_metadata.sh` fetches the ISA-Tab and the checksum lists and sizes every raw file without downloading it. `03_sample_table.sh` writes the MassCube sample tables (`config/sample_table_pos.csv`, `config/sample_table_neg.csv`) from the ISA-Tab alone, with `is_qc` and `is_blank` as "yes" and "no" strings, the only form MassCube reads; a column of booleans would silently be read as no QCs at all.

**One raw file at a time.** The deposit's 72,715,682,516 bytes are far more than the 10 GB this project allows itself, so `04_process_mode.sh` streams it: download one raw file with resumable parallel range requests, verify its sha256, convert it with `ThermoRawFileParser -f 1 -m 0 -z` (mzML, metadata as JSON, zlib off, vendor centroiding on), check the mzML (scan counts, every spectrum centroided and of the expected polarity, instrument, start and last scan time, and MassCube's own reading of its header), run MassCube's `process_single_file`, then delete the raw file and the mzML. What stays per file is MassCube's feature table and its `.mzpkl` scan cache, which gap filling reads in place of the mzML. On QC1 and QC11 the mzML was 161,869,144 and 167,426,249 bytes, about half the raw file, the scan cache 14,351,828 and 15,885,495 bytes and the feature table 28,405,508 and 30,184,789 bytes; conversion took 21 and 35 s and MassCube 70 and 109 s, at a peak resident memory of 914 MB (`results/conversion_manifest.tsv`). MassCube's own workflow function wants every mzML of a project at once, so a driver calls its functions in its order with the same `Params`.

MassCube orders injections for drift correction by each mzML's `startTimeStamp`, and its reader cannot parse the seven fractional digits ThermoRawFileParser writes (`2022-10-13T20:05:39.4620643Z`): it returns no time, and from its source its workflow would then keep the files in the order of the sample table, which it sorts with the QCs first, and put them all in one batch, without a warning. The driver parses the same attribute itself and keeps MassCube's empty reading in the manifest's `masscube_start_time` column.

**Pilot rules and alignment.** Before a full mode, six pilot files (the lowest and highest QC numbers and, by the seed 20261109, two HC, one INR and one IR) are processed; `05_pilot_rules.sh` applies the pre-registered m/z tolerance rule to them and projects the mode's time, memory and disk. `06_align_mode.sh` orders the files by acquisition time, flags conditioning QCs, applies the retention time and drift correction rules to the QCs, runs MassCube's alignment with gap filling, groups isotopes, adducts and in-source fragments, records per-injection metrics, and then deletes the mode's scan caches and gzips its feature tables so that the other mode fits.

**Annotation and normalisation.** `07_annotate.sh` runs only after both modes are aligned and their caches deleted. It downloads the newest MassCube MS/MS database version from Zenodo (on 9 November 2026, record 19121645 of 19 March 2026, 3,328,768,216 bytes, md5 db8c07244306af9e7b288ce68286eb7e) and checks its md5, loads it once and records what that costs, and takes the antiretroviral spectra out of the GNPS Drug Library (Zenodo record 17232042) for the 20 drugs in `config/antiretrovirals.tsv`. Annotation levels follow the pre-registration: level 2 is an identity search above 0.7 entropy similarity with at least 3 matched peaks and the precursor within 10 ppm; level 3 an analog or failing identity match; level 4 MS1 only; level 1 is not available without standards run on this instrument. Normalisation is MassCube's: LOWESS where the drift rule allows it, then probabilistic quotient normalisation.

**QC, filters, statistics and the published claims.** `08_qc_filters.sh` applies the failed-injection rule, tests whether group is confounded with injection order, logs each feature filter with its threshold and the number removed (detection in at least 50 per cent of one group, QC detection above 70 per cent, QC RSD at most 30 per cent, D-ratio at most 50 per cent; no blanks were deposited), keeps one representative per MassCube feature group, and computes the antiretroviral control. `09_statistics.R` fits the primary limma model and the pre-registered sensitivity analyses, `10_*` the exploratory pathway analysis, classifier and PERMANOVA, `11_replication.R` the verdict for each published claim and the original rule (Student t at p < 0.05 with OPLS-DA VIP > 1), and `12_figures.R` the figures.

## Deviations from the pre-registration

1. 9 November 2026. Section 4 says each file's start time is read by MassCube's `get_start_time` while its mzML exists. That function returns no time for ThermoRawFileParser's timestamps, so the driver parses the same `startTimeStamp` attribute, keeping six of its seven fractional digits, and MassCube's sorting and batch labelling are then applied to those times unchanged.
2. 12 November 2026. Section 4 orders feature grouping after annotation, as MassCube's workflow does. Grouping reads each feature's reference file from the scan caches, and silently leaves features ungrouped when a cache is missing; the caches cannot stay on disk until the database is downloaded (section 15). Grouping therefore runs at the end of alignment, before the caches are deleted, and a library match can no longer supply the adduct a feature is grouped as; MassCube's default ([M+H]+ or [M-H]-) is used instead.
3. 12 November 2026. Section 9 names `masscube.annotation.index_msp_to_pkl` for indexing the drug spectra. MassCube 1.2.22 has no such function; its `load_ms2_db` reads an MSP and indexes it in memory with the same unweighted entropy similarity it uses for the main database, and the antiretroviral subset is written as an MSP for it.
4. 12 November 2026. Section 13 has KEGG identifiers found through PubChem and the KEGG `conv` service. KEGG's `conv` maps PubChem substance identifiers, not compounds, so the identifiers are taken from the KEGG numbers PubChem lists among a compound's synonyms and confirmed in KEGG, with an exact KEGG name match as the fallback; the source of each is recorded in `config/kegg_compound_ids.tsv`.
5. 12 November 2026. The storage plan of section 15 did not foresee the size of MassCube's per-file feature tables, 28,405,508 and 30,184,789 bytes for the first two files and so, across a mode's 114 or 115 files, larger than its scan caches. They are gzipped at the end of alignment, after the per-injection metrics are taken from them.

Implementation choices the pre-registration left open, recorded here rather than as deviations: the PCA uses the features complete across the included injections; PERMANOVA needs complete data and uses half-minimum imputation for that analysis only; the classifier's inner loop is scikit-learn's `LogisticRegressionCV` over the pre-registered grid; and where ropls's automatic choice builds no OPLS-DA model, one orthogonal component is fixed and `results/oplsda_summary.tsv` says so.

## Repository structure

```
MassCube_HIV_metabolomics/
├── config/      sample tables, MassCube parameter files, conda specification and lock, recorded versions, published claims, antiretroviral list
├── scripts/     00 to 12 in order, plus lib/ with the shared Python, R and shell code
├── results/     every table the README cites (tracked)
├── figures/     every figure the README shows (tracked)
├── logs/        per-stage resource tables (tracked); console logs and stamps (not tracked)
├── data/        metadata, one raw file at a time, the databases; not tracked, see data/README.md
├── work/        MassCube project folders and transient mzML; not tracked
├── env/         the conda environment; not tracked
├── preregistration.md, run_all.sh, LICENSE, README.md
```

| script | does |
|---|---|
| `00_configure.sh` | detect or accept threads, RAM and disk, write `project.conf` |
| `01_install.sh` | conda environment, R packages, ThermoRawFileParser, versions recorded, caches emptied |
| `02_download_metadata.sh` | ISA-Tab, checksum lists and file sizes from MetaboLights |
| `03_sample_table.sh` | MassCube sample tables from the ISA-Tab |
| `04_process_mode.sh` | one raw file at a time: download, verify, convert, check, MassCube, delete |
| `05_pilot_rules.sh` | the pilot's m/z tolerance rule and the projection of the full mode |
| `06_align_mode.sh` | run order, retention time and drift rules, alignment, grouping, clean-up |
| `07_annotate.sh` | MassCube database and drug spectra, annotation levels, normalisation |
| `08_qc_filters.sh` | failed injections, run order test, filters, representatives, positive control, PCA |
| `09_statistics.R` | limma, sensitivity analyses, INR against IR, power |
| `10_permanova.R`, `lib/pathways.py`, `lib/classifier.py` | the exploratory analyses |
| `11_replication.R` | the published claims and the original rule |
| `12_figures.R` | the figures |

`scripts/lib/` also holds `common.sh` (configuration, stamps, the disk budget, resource logging), `measure.py` (peak resident memory without GNU time, which Windows lacks), `fetch.py` (resumable, checksummed downloads), `check_prose.sh`, `final_checks.sh`, `readme_numbers.R` (prints every number this README cites from its file), `check_dois.py` and `set_github_metadata.sh`.

## Usage

```bash
git clone https://github.com/Qasim-Hussain-Code/MassCube_HIV_metabolomics.git
cd MassCube_HIV_metabolomics
bash scripts/00_configure.sh --threads 16 --ram 16 --disk 10 --yes
bash run_all.sh
```

Or stage by stage, each of which skips when already done and says so:

```bash
bash scripts/01_install.sh
bash scripts/02_download_metadata.sh
bash scripts/03_sample_table.sh
bash scripts/04_process_mode.sh --mode pos --pilot
bash scripts/05_pilot_rules.sh --mode pos
bash scripts/04_process_mode.sh --mode pos
bash scripts/06_align_mode.sh --mode pos
# the same four for --mode neg, then:
bash scripts/07_annotate.sh
bash scripts/08_qc_filters.sh
bash scripts/lib/r_stage.sh 09_statistics scripts/09_statistics.R
bash scripts/lib/py_stage.sh 10_pathways scripts/lib/pathways.py
bash scripts/lib/py_stage.sh 10_classifier scripts/lib/classifier.py
bash scripts/lib/r_stage.sh 10_permanova scripts/10_permanova.R
bash scripts/lib/r_stage.sh 11_replication scripts/11_replication.R
bash scripts/lib/r_stage.sh 12_figures scripts/12_figures.R
```

Measured so far on the build machine (Windows 11 and Git Bash, configured by `00_configure.sh` with 16 threads and 16 GB of RAM):

| stage | elapsed | peak resident memory | note |
|---|---|---|---|
| 01 install | 2,007 s on the first build | not measured | the project grew from 0.00 to 2.43 GB (`logs/01_install_first_build.resources.tsv`) |
| 02 metadata | 55 s | 44 MB | `logs/02_download_metadata.resources.tsv`; its disk figures predate the removal of the doubly counted package cache |
| 03 sample tables | 2 s | 15 MB | `logs/03_sample_table.resources.tsv` |
| 04, per raw file | 21 and 35 s to convert, 70 and 109 s for MassCube | 914 MB | QC1 and QC11 (`results/conversion_manifest.tsv`); QC11's 350,562,367 bytes downloaded in 28 s |

The table for the full run, with every stage, replaces this one when it has run.

## Limitations

Those known before the results: the design is cross-sectional, with one centre; stool chemistry is shaped by diet and by ART, and every PLWH here was on ART; no per-sample covariate is public, so nothing is adjusted for; annotation is level 2 at best; no blanks were deposited, so contaminants and carry-over cannot be filtered against them; and the deposit's own metadata needed the corrections described under Data.

## Data availability

The inputs are public: MetaboLights MTBLS9264 (raw files, ISA-Tab and checksums), the MassCube MS/MS database (Zenodo record 19121645) and the GNPS Drug Library (Zenodo record 17232042). Nothing under `data/` or `work/` is tracked; `data/README.md` describes how the pipeline recreates it, file by file, against the published checksums.

## Citation

Method:

- Yu H, Ding J, Shen T, Liu M, Li Y, Fiehn O. MassCube improves accuracy for metabolomics data processing from raw files to phenotype classifiers. Nature Communications 2025, 16:5487. DOI 10.1038/s41467-025-60640-5, PMID 40593600.
- Broadhurst D, Goodacre R, Reinke SN, Kuligowski J, Wilson ID, Lewis MR, Dunn WB. Guidelines and considerations for the use of system suitability and quality control samples in mass spectrometry assays applied in untargeted clinical metabolomic studies. Metabolomics 2018, 14(6):72. DOI 10.1007/s11306-018-1367-3, PMID 29805336.

Data:

- Tian X, Xie Y, Yu L, Yao P, Dong M, Jin C, Wu N. Analysis of the gut microbiota and fecal metabolites in people living with HIV. Microbiology Spectrum 2024, 12(11):e00238-24. DOI 10.1128/spectrum.00238-24, PMID 39291988.
- Zhao HN, Kvitne KE, Brungs C, Mohan S, Charron-Lamoureux V, Bittremieux W, Tang R, Schmid R, et al. A resource to empirically establish drug exposure records directly from untargeted metabolomics data. Nature Communications 2025, 16:10600. DOI 10.1038/s41467-025-65993-5, PMID 41366203.

Software:

- Hulstaert N, Shofstahl J, Sachsenberg T, Walzer M, Barsnes H, Martens L, Perez-Riverol Y. ThermoRawFileParser: modular, scalable, and cross-platform RAW file conversion. Journal of Proteome Research 2020, 19(1):537-542. DOI 10.1021/acs.jproteome.9b00328, PMID 31755270.

Every DOI and PMID above was resolved through Crossref and Europe PMC on 12 November 2026 by `scripts/lib/check_dois.py` (`results/doi_check.tsv`). The full citation list (MetaboLights, ms_entropy and flash entropy search, limma, ropls, vegan, the annotation level schemes, PQN, Benjamini-Hochberg, Wieder et al. and the other method papers the pre-registration relies on) is added, and resolved the same way, with the results.

## License

The code is MIT (see `LICENSE`). MassCube is distributed by its authors under CC BY-NC 4.0, free for non-commercial use; it is installed from PyPI and not redistributed here. The MTBLS9264 data are distributed by MetaboLights under the EMBL-EBI Terms of Use, the MassCube database and the GNPS Drug Library by their authors under CC BY 4.0 on Zenodo, and the paper's supplement under CC BY.

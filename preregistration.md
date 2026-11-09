# Pre-registration: MassCube reanalysis of MTBLS9264

Written and committed on 9 November 2026, before any spectrum of MTBLS9264 was downloaded. It is not edited after it is committed and tagged `preregistration`. Anything decided later is recorded in the README under "Deviations from the pre-registration", with its date and reason.

## 1. What was known when this was written

The published results were known. Tian et al. (Microbiology Spectrum 2024, 12(11):e00238-24, doi 10.1128/spectrum.00238-24, PMID 39291988) and its supplement were read in full: 829 faecal metabolites, 292 called differential by two-tailed Student t-tests at p < 0.05, OPLS-DA in ropls on the two ion modes combined, four KEGG pathways (phenylalanine metabolism; arginine and proline metabolism; phenylalanine, tyrosine and tryptophan biosynthesis; arginine biosynthesis, at p = 0.005, 0.007, 0.0256 and 0.0457), and five named metabolites. N-acetyl-L-glutamate 5-semialdehyde, 11b-hydroxyandrost-4-ene-3,17-dione and 17alpha,21-dihydroxypregnenolone were reported higher in people living with HIV (PLWH); N-acetyl-L-phenylalanine and trehalose lower.

The spectra were not known. What had been examined is listed exactly:

- the ISA-Tab files `i_Investigation.txt`, `s_MTBLS9264.txt` and the two assay tables, each checked against the sha256 that MetaboLights publishes in `HASHES/metadata_sha256.json`;
- the FTP listing of `FILES/`, with the size of every raw file from a HEAD request and its published sha256 from `HASHES/data_sha256.json`;
- the study record from the MetaboLights web service.

The two metabolite assignment files (`m_..._v2_maf.tsv`, sha256 144a125a... positive and 8b078b74... negative) were downloaded with the metadata and checksum-verified, and have not been opened. They hold the authors' annotations and are opened only after the primary analysis below is complete, for the replication table of section 12.

## 2. Question and primary contrast

Do the faecal metabolomes of PLWH and healthy controls (HC) differ, feature by feature, when the deposited raw spectra are processed with MassCube and tested under the rules below, and which of the published claims survive that?

The primary contrast is PLWH against HC, within each ion mode separately. PLWH are the 32 immune non-responders (INR) and 38 immune responders (IR) of the deposit's disease factor, 70 in all; HC are 34. The INR against IR comparison is exploratory (section 13).

## 3. Data, as established from the metadata

| | Positive mode | Negative mode |
|---|---|---|
| Raw files on the FTP | 115 | 114 |
| Bytes | 38,488,887,662 | 34,226,794,854 |
| HC, INR, IR | 34, 32, 38 | 34, 32, 38 |
| Pooled QC files | 11 (QC1 to QC11) | 10 (nQC1 to nQC10) |
| Blanks | none deposited | none deposited |
| Converted mzML | none deposited | none deposited |

No per-sample covariate is deposited: no age, sex, ART regimen, CD4 count or viral load. Table S1 of the paper gives only group summaries. The models therefore carry no covariates.

Two source problems were found and decide how the tables are built.

1. In the positive assay table the Sample Name column is out of step with every other identifier column in all 114 rows. Extract Name, MS Assay Name, Derived Spectral Data File and Normalization Name all equal the raw file's name; Sample Name holds the same 114 names in another order, the ten QCs first, and read literally would put 40 files in the wrong group. The negative table agrees in every column. Files are therefore linked to participants through Extract Name, which agrees with the file names and their group prefixes in both modes. The check is rerun on every build and written to `results/isa_assay_consistency.tsv`. Section 9 adds an empirical test of the link.
2. `QC11.raw` is on the FTP with a published sha256 but in neither assay table, and the study table lists QC1 to QC10 only. It is included as a positive-mode pooled QC if its scans are all of positive polarity; otherwise it is excluded. A sensitivity analysis drops it (section 10).

The ISA protocols and the paper disagree on the instrument (Q Exactive against Q Exactive HF-X), resolving power (70,000 against 60,000), scans per cycle (10 against 8), column length (100 against 150 mm), flow (0.3 against 0.25 mL/min) and gradient (12 against 14 and 17 min). The raw files decide: instrument model, resolution where recorded, and the last scan time are read from each file and reported. The ISA extraction protocol adds 2-amino-3-(2-chlorophenyl)propionic acid at 4 ppm to every sample; it is used as a process internal standard (section 6).

## 4. Conversion and MassCube processing

Software: ThermoRawFileParser 1.4.5 (GitHub release, zip sha256 be1f8fd6f85b20750d5f59c324f031882f911035ff4e93c4e98c9b2d73944b37) and MassCube 1.2.22 (PyPI wheel sha256 ddb23d3b83d8b20da496cf3bb4ce3b5be94da39e1fcba51e2c928c60338a0ce6) with ms_entropy 1.3.4, in the environment recorded in `config/`.

Conversion: `ThermoRawFileParser -i <raw> -o <dir> -f 1 -m 0 -z`, that is, mzML, metadata as JSON, zlib compression off, and the vendor's peak picking left on, so spectra are centroided. After each conversion the file is checked: MS1 and MS2 scan counts, every spectrum flagged as centroid, every scan of the expected polarity, instrument model, run start time (`startTimeStamp`) and last scan time. A file that fails any check stops the run; it is not processed silently. Each file's raw sha256 is checked against the published value before conversion. MassCube's conversion guide recommends ProteoWizard's MSConvert with vendor peak picking and without zlib compression. ThermoRawFileParser is used instead because it runs unattended from a script on Windows and Linux alike; its default peak picking, like MSConvert's vendor filter, takes the centroids from Thermo's own library, and zlib is switched off as the guide asks.

MassCube parameters are those in `config/masscube_parameters_pos.csv` (sha256 f18dd1c6...) and `config/masscube_parameters_neg.csv` (06ea7f61...), which differ only in `ion_mode`. Peak height is the quantity; peak area is a sensitivity analysis. The intensity thresholds are the values MassCube's parameter definitions recommend for Orbitrap data, 30,000 for MS1 and 10,000 for MS2 (its automatic default for an Orbitrap file, used only when no parameter file is given, sets 50,000 for MS1); the m/z and retention time tolerances are MassCube's defaults (0.01 Da for MS1 and alignment, 0.015 Da for MS2, 0.2 min for alignment). Two values in those files start as `false` and may be switched on by the rules of section 7: `correct_rt` and `signal_normalization`. Fuzzy (analog) search is on, the MS2 similarity threshold is 0.7, `plot_normalization` is off, and MassCube's own `run_statistics` (a UMAP plot) is off.

The storage plan of section 15 means that no mode's mzML files exist together, while MassCube's `untargeted_metabolomics_workflow` expects every mzML of a project in its data folder at once. The workflow function is therefore not called as a whole. A driver calls the same functions the workflow calls, in the same order and with the same `Params` object: `process_single_file` for each file as soon as it is converted (writing the per-file feature table and the `.mzpkl` scan cache that gap filling reads), then `feature_alignment` with gap filling, the pickle of aligned features, `annotate_aligned_features` once the database is present, `group_features_after_alignment`, the table and MSP exports, `signal_normalization` when switched on, then `sample_normalization`. The one replaced step is reading each file's start time: MassCube's own `get_start_time` is called on each mzML while it exists and the value is stored, instead of being read again later. Injection order is the order of those start times, and the driver checks that it equals MassCube's `analytical_order`.

## 5. Injection order and confounding

Injection order is each file's `startTimeStamp`, ties broken by file name. Batches are MassCube's `label_batch_id` assignment and are reported.

Whether group is confounded with order is tested in each mode with a two-sided Wilcoxon rank-sum test of injection order between PLWH and HC samples, QCs excluded. A p value below 0.05 is reported as order associated with group, and the sensitivity analysis of section 10 that adds injection order to the model then becomes part of the headline. PC1 and PC2 of the QC-included PCA (section 8) are tested against order by Spearman correlation, study samples only.

QC injections that come before the first study sample, other than the last of them, are treated as conditioning injections (Broadhurst et al. 2018): they are kept in the tables, flagged, and left out of RSD, D-ratio and drift correction.

## 6. Failed injections

An injection is excluded from the primary analysis if any of the following holds, each computed within its mode:

- MassCube's own check flags it: its highest feature is below one third of the median of that value across files (`run_evaluation`);
- its feature count, or the sum of its feature heights, is below one third of the mode median;
- its last MS1 scan time is below 90 per cent of the mode median, that is, the run was truncated;
- the internal standard (positive [M+H]+ m/z 200.0473, negative [M-H]- m/z 198.0327, within 10 ppm) is not detected. This criterion is dropped for a mode if the standard is found in fewer than 90 per cent of that mode's files, since it would then measure the standard rather than the injection.

Every injection's values and its verdict are written to `results/injection_qc.tsv`. The analysis is repeated with excluded injections put back (section 10).

## 7. Rules decided by the QC data before alignment

All three use QC data only and are applied before any study sample is compared with another.

m/z tolerance, decided on the pilot of section 16. The 200 most intense features of the pilot file of median total intensity are matched in every other pilot file within 0.02 Da and 0.5 min. If the 99th percentile of the per-feature m/z range exceeds 0.01 Da, `mz_tol_ms1` and `mz_tol_alignment` are set to 0.015 Da for that mode and the pilot is reprocessed. The pilot's six raw files are kept until this rule has been applied, so that reprocessing needs no second download, and are deleted afterwards. Tolerances are never narrowed. The internal standard's mass error in ppm is reported either way.

Retention time correction, decided after single-file processing of a whole mode. The same matching is done across the non-conditioning QCs. `correct_rt` is switched on if the median per-feature retention time range across QCs exceeds 0.05 min.

Drift correction. MassCube's QC-based LOWESS is applied in a mode only if, after conditioning injections are set aside, at least 8 QCs remain (Broadhurst et al. 2018), a QC comes before the first study sample and another after the last, and no run of more than 20 consecutive study injections lacks a QC. Otherwise the mode is left uncorrected and the README says so. If no QC remains in a mode at all, RSD, D-ratio and LOWESS are dropped for it and that goes first in Limitations.

## 8. Normalisation, transformation and filters

Order as MassCube applies it: LOWESS drift correction where section 7 switches it on, then probabilistic quotient normalisation through MassCube (`sample_norm_method = pqn`), then log2. Values that are still zero after gap filling become missing. None is imputed for the primary test; half of the feature's minimum observed value is the imputation for a sensitivity analysis.

Feature filters, applied in this order, each logged with its threshold, the number removed and the reason in `results/<mode>_filter_log.tsv`:

1. Detected before gap filling (`feature_id_arr` not -1 in MassCube's aligned features) in at least 50 per cent of the samples of at least one group, HC or PLWH.
2. Detected before gap filling in more than 70 per cent of the non-conditioning QCs (Broadhurst et al. 2018).
3. QC relative standard deviation at most 30 per cent, on normalised peak heights. Broadhurst et al. give 20 per cent for some sample types and 30 per cent for others; stool is among the most heterogeneous matrices, so 30 is the primary threshold and 20 a sensitivity analysis.
4. D-ratio at most 50 per cent: the standard deviation of log2 QC values over that of log2 study-sample values.
5. Blank filter: not applicable, no blanks were deposited.

Isotopes, adducts and in-source fragments that MassCube places in one `feature_group_id` are then collapsed to one representative for testing: among the group's features that passed the filters and are not flagged as isotopes or in-source fragments, the one with the highest median normalised height across study samples. Features with no group are their own representative.

PCA is run on log2, autoscaled data with QCs included and drawn coloured by group and by injection order. Tight QCs are the check.

## 9. Annotation and the antiretroviral positive control

The MassCube MS/MS database is loaded once, after both modes are processed and their spectra deleted, and both modes are annotated from their saved aligned features. The version is the newest on Zenodo at download time; on 9 November 2026 that is record 19121645 (19 March 2026, `MassCube_DB_2026-03-19.pkl`, 3,328,768,216 bytes, md5 db8c07244306af9e7b288ce68286eb7e), which was built with ms_entropy 1.3.4. The pickle must load before it is relied on. Search mode, similarity and matched peak number are kept for every hit.

Confidence levels, after Sumner et al. 2007 and Schymanski et al. 2014:

- Level 2: MassCube identity search (`search_mode` = `identity_search`) with entropy similarity above 0.7 (MassCube reports a hit only when the similarity exceeds `ms2_sim_tol`), at least 3 matched peaks, and a precursor m/z within 10 ppm of the matched library precursor.
- Level 3: any fuzzy (analog) search hit, any identity hit failing a level 2 criterion, or a compound class assignment.
- Level 4: MS1 only.
- Level 1 needs authentic standards run on this instrument and is not available.

Only level 2 annotations are named as compounds in the README; level 3 is reported as such.

The positive control is antiretroviral drugs. All 70 PLWH had been on ART for two years with undetectable viral load (ISA protocol), so drug-derived features should appear in PLWH and not in HC. The drugs looked for are tenofovir and its prodrugs tenofovir disoproxil and tenofovir alafenamide, emtricitabine, lamivudine, efavirenz, dolutegravir, bictegravir, raltegravir, elvitegravir, rilpivirine, nevirapine, lopinavir, ritonavir, darunavir, cobicistat, zidovudine and abacavir, with azvudine and ainuovirine, both approved in China. Where the MassCube database lacks a reference spectrum, the GNPS Drug Library (Nature Communications 2025, doi 10.1038/s41467-025-65993-5) is used, indexed for only those compounds with `masscube.annotation.index_msp_to_pkl`. A drug counts as detected in a sample when a level 2 feature for it was detected there before gap filling.

The control passes if at least one antiretroviral is detected in at least 50 per cent of PLWH and in at most one HC sample. Detection by group is reported for every drug, and any HC detection is reported as a possible label, carry-over or annotation problem. The same counts are computed under the Sample Name link that section 3 rejected, as an empirical test of the link: the link under which drugs are confined to PLWH is the one the biology supports, and if that is not the Extract Name link, the analysis stops for review.

Drug-derived features (parent drugs, and features whose level 2 or 3 match is to a drug or a known drug metabolite) are a class of their own. They are reported, kept out of the endogenous pathway analysis and out of the main classifier. Each of the five published named metabolites is checked against the drug features for an explanation by a drug: formula, adduct, and MS2 similarity to the drug spectra.

## 10. Primary test and sensitivity analyses

Primary test, per mode, in limma on log2 normalised heights of the representative features: design `~ group` with HC as reference, `lmFit`, `eBayes` with its defaults (no intensity trend, not robust), two-sided, Benjamini-Hochberg within the mode, significant at q < 0.05. A feature is tested only if it has at least 3 non-missing values in each group. The effect is the log2 fold change, PLWH over HC, with its 95 per cent confidence interval from `topTable(confint = 0.95)`. No covariates, since none are public per sample.

Sensitivity analyses, each reported as the number of primary hits it keeps and the sign agreement:

- Welch t-test, and Wilcoxon rank-sum test, each with Benjamini-Hochberg within mode;
- with and without LOWESS drift correction (whichever the primary did not use);
- with excluded injections put back;
- peak area instead of peak height;
- half-minimum imputation of missing values;
- QC RSD at most 20 per cent instead of 30;
- injection order (linear and quadratic terms) added to the limma design;
- positive mode without QC11.

## 11. Power

The minimum standardised difference detectable at 80 per cent power, two-sided, for the realised group sizes, by the exact power of a two-sample t-test, computed by integrating the normal tail probabilities over the chi-square distribution of the variance. statsmodels' `TTestIndPower` was tried first and is not used: on this build it gave a smaller detectable difference at alpha 0.05/1,000 than at 0.05/500, and scipy's noncentral t, which it relies on, returned undefined values at the noncentralities involved. For 70 against 34: d = 0.591 at alpha 0.05; at alpha 0.05 divided by the number of filtered features, d = 1.027 for 500 features, 1.067 for 1,000, 1.105 for 2,000 and 1.154 for 5,000. The value for the realised number of filtered features and the realised group sizes after exclusions is computed by the statistics stage and written to `results/power.tsv`.

## 12. The published claims

`results/replication_tian2024.tsv` has one row per claim: the five named metabolites, the four pathways, and the headline counts 829 and 292. Each row gives whether the claim was detected here, best annotation level and mode, direction, log2 fold change with its interval, p and q, a verdict and a one-line reason.

Named metabolite, matched by InChIKey (first block if stereochemistry differs) to a level 2 or level 3 annotation in either mode; where several features match, the one with the best level, then the lowest q. The InChIKey of each published name is taken from PubChem by name and recorded with its CID in `config/tian2024_claims.tsv` before the authors' assignment files are opened. Isomers that share a formula and cannot be told apart by MS2 (trehalose among the disaccharides) are reported as such rather than counted as a match.

- held: level 2, q < 0.05 in the primary test, same direction as published;
- partly held: same direction and either level 3 with q < 0.05, or level 2 or 3 with p < 0.05 and q of 0.05 or more;
- not supported: annotated at level 2 or 3, and either p of 0.05 or more or the opposite direction;
- not testable: not annotated at level 3 or better, or removed by the filters of section 8.

Pathway, by the over-representation test of section 13 against KEGG human pathways:

- held: q < 0.05;
- partly held: p < 0.05 and q of 0.05 or more;
- not supported: p of 0.05 or more with at least 3 background members;
- not testable: fewer than 3 background members.

Headline counts, compared with the counts this reanalysis gives under the original rule below (829 against the number of level 2 and 3 annotated compounds after filtering, de-duplicated by InChIKey first block across modes; 292 against the number the original rule calls among them): held if within 25 per cent of the published count, partly held if within 50 per cent, otherwise not supported.

The original rule is rerun on these data so that a reader can see how much of any difference comes from processing and how much from statistics: Student t-test with equal variances, unadjusted p < 0.05, and OPLS-DA VIP > 1, both modes combined, on the annotated compounds. OPLS-DA in ropls with one predictive component, unit-variance scaling, 7-fold cross-validation and 1,000 permutations; R2Y, Q2 and both permutation p values are reported, and VIP is described as descriptive, not inferential.

## 13. Exploratory analyses

Labelled exploratory wherever they appear.

- PERMANOVA: vegan `adonis2`, Euclidean distance on autoscaled log2 data, 9,999 permutations, per mode, for comparison with the paper's 10,000.
- Classifier, exploratory and not diagnostic: elastic-net logistic regression in scikit-learn on the representative features of both modes, drug-derived features excluded, under repeated nested cross-validation (5 outer folds, 5 inner, 10 repeats, stratified). Inside each training fold: half-minimum imputation, log2, autoscaling; the inner loop chooses `l1_ratio` from 0.1, 0.5 and 0.9 and `C` from 12 values log-spaced from 0.001 to 100 by AUC. Each sample's out-of-fold probability is averaged over the 10 repeats; AUC is computed on those, with a 95 per cent interval from 2,000 stratified bootstrap resamples. Negative control: the same procedure with permuted labels. The same model with drug-derived features included shows how much ART exposure alone separates the groups.
- Pathway over-representation, after Wieder et al. 2021: Fisher's exact test, one-sided for over-representation, of endogenous level 2 compounds significant in the primary test (either mode, de-duplicated by InChIKey first block), against a background of all endogenous level 2 compounds detected after filtering, never all of KEGG. KEGG compound identifiers come from the MassCube database record where present, otherwise from InChIKey through PubChem and the KEGG REST `conv` service; unmapped compounds are counted and reported. KEGG human (hsa) pathways for comparison with the paper; KEGG reference pathways as a sensitivity analysis, since stool chemistry is partly microbial. Benjamini-Hochberg across pathways with at least 3 background members. Membership is fetched from the KEGG REST API and cached in `config/` with the date.
- INR against IR, with the limma model of section 10 restricted to PLWH.

## 14. Seed

20261109 wherever randomness enters: pilot selection, permutations, cross-validation folds, bootstrap and figure jitter.

## 15. Storage plan and its fallbacks, in order

The project, environment and caches included, must stay under 10 GB, and at least 2 GB must stay free on the drive (C: on the build machine, which Windows shares with its pagefile). Before every download the project size plus the incoming file must stay under 9.5 GB, and the free space minus the incoming file above 2 GB; otherwise the run stops with a plain message.

1. Primary plan, one ion mode at a time, one file at a time: download the raw file, verify its sha256, convert it, check the mzML, run MassCube's single-file processing, record its start time, delete the raw file and the mzML. MassCube's per-file feature table and `.mzpkl` scan cache are kept, since gap filling reads the cache and not the mzML. Then align the mode, save its aligned features, and delete its `.mzpkl` files; the manifest, checksums and timestamps that regenerate them stay. Then the other mode. Only after both are aligned and their caches deleted is the 3.3 GB database downloaded, and both modes are annotated from their saved aligned features.
2. If the pilot projects that one mode's `.mzpkl` caches would not fit, each cache is rewritten after single-file processing with the peak arrays of its MS2 scans emptied. Gap filling reads only MS1 scans from the cache, so its values should not change; this is checked on the pilot files by running gap filling with full and reduced caches and comparing every filled value, and the switch is made only if they are identical. Compressing the mzML with zlib is not among the fallbacks: only one mzML exists at a time, so its size does not decide whether a mode fits, and MassCube's conversion guide says not to use zlib compression.
3. If that still cannot fit, the run stops and the author is asked.
4. The database is never on disk together with a mode's scan caches.

## 16. Pilot

Before a full mode, six files of that mode are processed through the whole per-file path: two QCs (the lowest and highest QC numbers) and, chosen with the seed, two HC, one INR and one IR. Time, peak resident memory and disk per file are measured and the full mode is projected; the projection and later the actual figure are recorded. Files are processed one at a time in a single process as they are converted, so MassCube's own worker count does not apply to single-file processing. Peak resident memory of every step must stay under 14 GB on a 16 GB machine and below the memory actually free on the build machine, where other applications hold much of it; if the pilot or alignment exceeds that, the run stops for review.

## 17. What will be reported regardless of outcome

Every item above, with its rule applied as written: the positive control whether it passes or fails, every replication verdict including the unfavourable ones, every sensitivity analysis, and every failed or excluded file. Nothing is dropped for giving an inconvenient answer. If a rule turns out to have been badly chosen, its result is reported first and the reasoning about it follows, labelled as a deviation with its date.

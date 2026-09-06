# Reproducibility protocol

## Reproducibility target

The target is a complete, deterministic path from the supplied aggregate annual
edge lists and node registry to:

- validated empirical and known-truth simulation tables;
- fitted-kernel and structural-role objects;
- all publication figures, including the graphical abstract;
- the HTML review manuscript and MethodsX Word manuscript;
- a one-page MethodsX cover letter;
- session, configuration, input, figure, and output provenance records.

The publication configuration uses seed `20260817`, 200 simulation replicates,
1,000 empirical random cut-floor samples, 128 simulation cut-floor samples, and
200 alternating-maximization restarts per sign.

## Execution levels

### 1. Diagnostic gate

Use this after code or data changes. It performs a clean analysis with smaller
simulation and cut-search settings and does not render the manuscript.

```r
Sys.setenv(
  METHODSX_MODE = "diagnostic",
  METHODSX_REUSE_INTERMEDIATE = "0",
  METHODSX_RENDER_RMD = "0",
  METHODSX_SKIP_SIMULATION = "0",
  METHODSX_CUT_RESTARTS = "20"
)
source("run_all.R")
```

### 2. Full publication reproduction

Use this for a release candidate and after changes affecting analysis code,
inputs, simulation design, or estimator settings.

```r
Sys.setenv(
  METHODSX_MODE = "publication",
  METHODSX_REUSE_INTERMEDIATE = "0",
  METHODSX_RENDER_RMD = "1",
  METHODSX_SKIP_SIMULATION = "0",
  METHODSX_CUT_RESTARTS = "200"
)
source("run_all.R")
```

### 3. Manuscript-only publication render

Use this after prose, title, author metadata, caption, or graphical-abstract
layout changes when publication intermediates have already passed validation.

```r
Sys.setenv(
  METHODSX_MODE = "publication",
  METHODSX_REUSE_INTERMEDIATE = "1",
  METHODSX_RENDER_RMD = "1",
  METHODSX_SKIP_SIMULATION = "1",
  METHODSX_CUT_RESTARTS = "200"
)
source("run_all.R")
```

## Required acceptance checks

1. `METHODSX_PIPELINE_validation.csv` contains 16 rows and every `passed` value
   is `TRUE`.
2. The maximum calibrated null false-positive rate is no greater than 0.10; the
   current publication result is 0.05.
3. The empirical node universe contains exactly 175 aligned nodes for every
   year from 2006 through 2024.
4. The 2016 and 2019 total weights are exactly 358367 and 384035.
5. All 18 consecutive transitions have finite evidence scores.
6. All three estimators and all four fixed resolutions are represented.
7. The figure manifest contains 16 current PDF/PNG files.
8. The output manifest contains the current HTML and DOCX products and excludes
   its own path.
9. The manuscript DOCX contains no updateable fields or external file targets;
   ordinary DOI and e-mail hyperlinks may remain.
10. `output/Cover_Letter_MethodsX.docx` exists and its declarations have been
    confirmed by both authors before submission.
11. `scripts/release_check.R` completes without error.
12. A human reviewer inspects the complete HTML and DOCX, with special attention
    to figure widths, the retained template header, mathematical notation, and
    the graphical abstract at thumbnail scale.

## Cut-norm terminology

The implementation uses alternating maximization on both signs, deterministic
and random restarts, and a legacy random subset-pair floor. It reports a
reproducible lower-bound estimate, not an exact cut norm. Exact cut-norm
calculation is NP-hard, and the reported node-aligned quantity is not an
unlabeled graphon cut distance because no relabeling optimization is performed.

## Expected computational claims

The manuscript's numerical statements are tied to the publication summaries.
Key fixed checkpoints include:

- maximum null false-positive rate: 5%;
- medium-effect raw cut lower-bound detection: 86.0% for directional shifts,
  61.5% for local block shocks, and 59.5% for role reassignments;
- strongest empirical six-perspective consensus: 2017–2018, 84.8%, with five
  of six perspectives at or above the 75th percentile;
- 2015–2016 ranks eleventh by consensus and is not a dominant break.

If any checkpoint changes after a clean publication run, revise the narrative
from the new CSV outputs before release.

## Environment capture

`DESCRIPTION` declares the direct package dependencies and
`METHODSX_PIPELINE_sessionInfo.txt` records the exact environment used for the
publication run. The GitHub Actions workflow performs a clean publication run
on a second operating system. Cross-platform byte-identical graphics are not
guaranteed because font and graphics devices may differ; analytical CSV values
and validation outcomes are the primary cross-platform invariants.

## Word field and cover-letter handling

`scripts/production/07_freeze_docx_fields.R` first creates
`output/MethodsX-reference-clean.docx` from the unchanged
`MethodsX-reference.docx`. The runtime copy replaces the template's named
`rIdMethodsXHeader` relationship with the next available numeric relationship
ID and updates the matching header reference. This prevents the repeated
`header1.xml` relationship/coercion warnings from `officer` and `officedown`
without removing or rebuilding the MethodsX header.

After rendering, the same script converts the seven officedown figure-number
`SEQ` fields to fixed integers. This removes Word's generic external-field
update prompt without unlinking embedded images or removing DOI/e-mail
hyperlinks. DOCX archives are extracted through `zip::unzip()` from normalized
absolute paths and are checked for `word/document.xml` before processing. The
release check reopens both the clean reference and final manuscript and rejects
non-numeric relationship IDs, remaining field markers, update-on-open settings,
or external file targets.

`scripts/production/08_cover_letter.R` creates
`output/Cover_Letter_MethodsX.docx`. Its date defaults to the run date; for a
specific submission date use, for example,
`Sys.setenv(METHODSX_SUBMISSION_DATE = "2026-08-28")`. Confirm the author-
approval, exclusivity, originality, and conflict-of-interest statements before
the letter is sent. The release ZIP also includes the visually reviewed copy at
`submission/Cover_Letter_MethodsX.docx`.

## Release sequence

1. Run the diagnostic gate.
2. Run the full publication reproduction.
3. Run `source("scripts/release_check.R")`.
4. Inspect HTML and DOCX visually.
5. Inspect and approve the generated cover letter.
6. Have the second author validate software and computational claims.
7. Confirm the CRediT statement and final author approval.
8. Select code and data licenses.
9. Tag the release and attach the reproduction ZIP plus the graphical abstract.

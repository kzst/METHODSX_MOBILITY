# METHODSX_MOBILITY

Reproducible materials for **Comparing Network Indicators and Graphon Distances
for Structural Anomaly Detection: A Reproducible Workflow for Directed Weighted
Longitudinal Networks**.

Authors: Zsolt Tibor Kosztyán
([ORCID 0000-0001-7345-8336](https://orcid.org/0000-0001-7345-8336)) and
Kornél Dénes
([ORCID 0009-0002-9527-6124](https://orcid.org/0009-0002-9527-6124)).

## What this repository does

The workflow compares and integrates two complementary sources of evidence for
change in node-aligned directed weighted longitudinal networks:

1. conventional annual network indicators describing activity, connectivity,
   reciprocity, and concentration; and
2. graphon-inspired fitted-kernel distances, adjacency-spectrum benchmarks, and
   structural-role instability.

Known-truth simulations calibrate channel-specific detection thresholds. The
empirical layer converts consecutive-year changes to within-channel anomaly
scores, summarizes correlated channels as six transparent perspectives, and
reports convergence and disagreement rather than forcing a single metric to
dominate.

## Reproduced application

The included aggregate annual networks represent Hungarian higher-education
applications among 175 micro-regions from 2006 through 2024. Edges are annual
aggregate application counts, not individual-level trajectories and not
realized enrolments.

## Quick start

Requirements:

- R 4.4.0 or newer;
- Pandoc (bundled with recent RStudio releases);
- the R packages declared in `DESCRIPTION`.

Install missing packages:

```r
source("scripts/install_dependencies.R")
```

Run a fast diagnostic analysis without manuscript rendering:

```r
Sys.setenv(
  METHODSX_MODE = "diagnostic",
  METHODSX_REUSE_INTERMEDIATE = "0",
  METHODSX_RENDER_RMD = "0",
  METHODSX_CUT_RESTARTS = "20"
)
source("run_all.R")
```

Run the complete publication pipeline from raw aggregate inputs:

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

Re-render the manuscript from validated publication intermediates:

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

## Expected checks and products

A successful publication run reports all 16 analytical validation checks as
`TRUE` and verifies 16 figure files: PDF and PNG exports for the graphical
abstract, Figures 1–6, and Supplementary Figure S1.

Principal outputs:

- `output/MethodsX_MOB.docx`
- `output/Cover_Letter_MethodsX.docx`
- `output/MethodsX_MOB_figures.html`
- `output/figures/FIG00_graphical_abstract.pdf`
- `output/figures/FIG00_graphical_abstract.png`
- `output/diagnostics/METHODSX_PIPELINE_validation.csv`
- `output/diagnostics/METHODSX_TEMPLATE_compliance.csv`
- `output/diagnostics/METHODSX_PIPELINE_output_manifest.csv`
- `output/diagnostics/METHODSX_FIGURE_manifest.csv`

The output manifests record byte sizes and MD5 hashes after rendering. The
output manifest deliberately excludes itself. Word rendering uses the exact
MethodsX Method Article Template supplied through the editor's link. After
rendering, automatic caption fields become fixed text, template comments are
removed, all article images remain embedded, DOI/e-mail hyperlinks remain
clickable, and Word does not need to ask whether fields referring to other
files should be updated. A separate 19-check structural audit verifies the
mandatory template headings and order, title and abstract limits, completed
Specifications table, author ORCIDs, separate graphical abstract, declarations,
and removal of template instructions.

## Repository map

```text
METHODSX_MOBILITY/
├── MethodsX_MOB.Rmd
├── MethodsX-Method-Article-Template.docx  # exact editor-linked template
├── run_all.R
├── DESCRIPTION
├── CITATION.cff
├── submission/
│   ├── Cover_Letter_MethodsX.md    # editable text and declaration gate
│   ├── AUTHOR_METADATA.csv          # Editorial Manager author-data sheet
│   ├── EDITORIAL_MANAGER_CHECKLIST.md
│   └── RESPONSE_TO_EDITOR.md
├── data/
│   ├── raw/                 # aggregate annual inputs and spatial reference
│   └── derived/             # deterministic analytical products
├── scripts/
│   ├── production/          # active pipeline and graphical-abstract code
│   ├── reference/revised/   # graphon-distance reference implementation
│   ├── install_dependencies.R
│   └── release_check.R
└── output/                  # manuscripts, figures, diagnostics, and objects
```

Detailed execution and provenance guidance is in `REPRODUCIBILITY.md`; field
definitions are in `DATA_DICTIONARY.md`; the conceptual-figure storyboard is in
`GRAPHICAL_ABSTRACT_DESIGN.md`. The cover letter is generated reproducibly by
`scripts/production/08_cover_letter.R`; set `METHODSX_SUBMISSION_DATE` to an
ISO date (`YYYY-MM-DD`) before rendering if the submission date should differ
from the current date.

## Authorship and validation gate

Kornél Dénes is assigned the CRediT roles Software, Validation, and Writing –
review and editing. Before submission, the independent code validation must be
completed, both authors must approve the final manuscript, and both must accept
accountability for their contributions. Editorial Manager must contain an
institutional e-mail address and the exact ORCID for each author, and both
authors must complete the platform's authorship verification. The institutional
addresses are `kosztyan.zsolt@gtk.uni-pannon.hu` (Zsolt Tibor Kosztyán) and
`morgosz@student.elte.hu` (Kornél Dénes); both are recorded in
`submission/AUTHOR_METADATA.csv`. The validation protocol is
recorded in `CONTRIBUTING.md`, and the resubmission gate is recorded in
`submission/EDITORIAL_MANAGER_CHECKLIST.md`.

## Data responsibility

Only aggregate micro-region-level edge counts and node metadata belong in the
public reproducibility package. Do not add applicant-level records, direct
identifiers, credentials, local Dropbox metadata, or restricted administrative
files to the repository.

## License before public release

No software or data license is imposed by this preparation step. The repository
owner must select compatible code and data licenses before making the GitHub
repository public. In the absence of a license, reuse rights are not granted by
default.

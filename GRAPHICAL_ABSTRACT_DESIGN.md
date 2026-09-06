# Graphical abstract design

## Purpose

The graphical abstract must make one claim visible at a glance: **graphon-inspired
distances complement, rather than replace, conventional longitudinal network
analysis**. Agreement between the two evidence families supports a broader
structural-anomaly claim; disagreement reveals the type, scale, or specification
sensitivity of the change.

## Evidence-based design principles

The design follows the guidance in Chapters 7, 8, and 10 of the author's
`RESBOOK_EN` manuscript:

- one image and one take-home message;
- minimal text that remains legible at thumbnail size;
- a left-to-right input → comparison → validation → contribution flow;
- an explicit link between research questions, methods, and contribution;
- a consistent colorblind-safe palette;
- accurate diagram-as-code rather than decorative generative imagery;
- human verification of every label and relationship.

## Storyboard

| Stage | Visual content | Scientific role |
|---|---|---|
| Longitudinal networks | Annual directed weighted snapshots, 2006–2024, 175 aligned nodes | Defines the data and asks when, what type, and how robust |
| Annual network indicators | Activity, connectivity, reciprocity, concentration | Preserves familiar, directly interpretable yearly trajectories |
| Graphon-inspired distances | Raw intensity, unit-mean shape, Frobenius, cut lower bound, spectral modes, structural roles | Retains coordinated relational information that scalar summaries may miss |
| Calibration and robustness | Known-truth simulations, within-channel anomalies, estimator/resolution sensitivity, six-perspective profile | Makes heterogeneous evidence comparable without hiding disagreement |
| Reproducible contribution | Convergence strengthens claims; disagreement identifies change type and sensitivity | States the method-selection and interpretation contribution |

Each stage also uses a compact, literal symbol that remains secondary to the
scientific text: a small directed network for the input, three indicator time
series with one red anomaly marker, a fitted graphon surface, a calibration
target with a check mark, and converging/diverging arrows for interpretation.
The symbols are drawn by the same R source as the boxes and labels; no external
image file is linked into the Word manuscript.

## Production choice

The graphical abstract is implemented in
`scripts/production/06_graphical_abstract.R` and called from
`MethodsX_MOB.Rmd`. The pipeline exports `FIG00_graphical_abstract.pdf` and
`FIG00_graphical_abstract.png`. This option was selected over a manually edited
design because it gives a single auditable source, deterministic regeneration,
vector output, consistent colors, and automatic inclusion in the figure and
MD5 manifests.

## Accessibility and layout constraints

- Landscape aspect ratio: 2:1.
- Word display width: 6.8 inches, kept within the template's text area.
- Primary colors: blue for conventional indicators, orange for graphon-inspired
  distances, gold for calibration, and olive for contribution.
- No color alone carries meaning; every lane and outcome is labeled.
- The graphical abstract is unnumbered so the article's Figure 1–6 sequence is
  unchanged.
- The three research questions in the input panel are split across two lines
  (`WHEN? WHAT?` / `HOW ROBUST?`) and held above the lower panel boundary so
  they remain visible at the 6.8-inch Word display size.
- The red anomaly guide is a bounded segment inside the indicator mini-chart;
  it must never extend through neighbouring panels or the footer.
- Icon detail is deliberately limited; if a future edit reduces text clearance,
  label visibility takes priority over decorative detail.

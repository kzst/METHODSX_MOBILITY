# Data dictionary

## Raw aggregate inputs

### `data/raw/empirical/nodes.txt`

| Field | Meaning |
|---|---|
| `nodeID` | Stable micro-region identifier used for node alignment |
| `nodeLabel` | Human-readable micro-region label |
| `nodeLat` | Latitude used for reference and visualization |
| `nodeLong` | Longitude used for reference and visualization |

### `data/raw/empirical/YYYY_from_to_kist_weight.txt`

| Field | Meaning |
|---|---|
| `from` | Applicant-residence micro-region identifier |
| `to` | Chosen-campus micro-region identifier |
| `weight` | Aggregate number of application records for the directed pair and year |

There is one file for each year 2006–2024. These are aggregate annual edge
counts. They do not contain applicant identifiers or individual pathways.

### `data/raw/spatial/hungary.RData`

Spatial reference object used by the retained reference workflow. The active
MethodsX figures do not require applicant-level or individual geolocation data.

## Derived empirical outputs

| File | Main content |
|---|---|
| `EMPIRICAL_annual_input_checks.csv` | Node counts, edge counts, total weights, and weight validity by year |
| `EMPIRICAL_annual_indicators.csv` | Eight interpretable network indicators for 19 annual snapshots |
| `EMPIRICAL_transition_evidence.csv` | Channel scores and anomaly evidence for 18 consecutive transitions |
| `EMPIRICAL_estimator_agreement.csv` | Rank and top-five agreement across fitted-kernel estimators |
| `EMPIRICAL_resolution_agreement.csv` | Rank and top-five agreement across fixed resolutions |
| `EMPIRICAL_role_memberships_K13.csv` | Fixed-resolution structural-role membership and role descriptors |
| `EMPIRICAL_role_transition_summary.csv` | Consecutive-year role instability summaries |
| `EMPIRICAL_multiresolution_refinement.csv` | Approximate refinement consistency across K = 2, 4, 8, 13 |
| `EMPIRICAL_role_flow_matrix_2018.csv` | Observed application-weight flows between destination-role ranks |
| `EMPIRICAL_destination_core_definition.csv` | Destination-core membership and reference fields |
| `EMPIRICAL_node_aligned_networks.rds` | Prepared node registry, graph objects, years, sources, and checks |

Structural roles group nodes with similar directed connectivity profiles. They
are not modularity communities and do not rank institutional quality.

## Derived simulation outputs

| File | Main content |
|---|---|
| `SIMULATION_scenario_specification.csv` | Null and six known perturbation designs with effect sizes |
| `SIMULATION_kernel_design.csv` | Expected before/after kernels for visual inspection |
| `SIMULATION_raw_scores.csv` | Replicate-level channel scores |
| `SIMULATION_null_thresholds.csv` | Channel-specific thresholds calibrated from the null |
| `SIMULATION_detection_summary.csv` | Calibrated detection rates by scenario, effect, and channel |

## Diagnostic and provenance outputs

| File | Purpose |
|---|---|
| `EMPIRICAL_input_manifest.csv` | Size and MD5 hash of every raw empirical input |
| `METHODSX_PIPELINE_configuration.csv` | Effective runtime configuration |
| `METHODSX_PIPELINE_validation.csv` | Deterministic acceptance checks |
| `METHODSX_PIPELINE_sessionInfo.txt` | R, platform, and package versions |
| `METHODSX_FIGURE_manifest.csv` | Size and MD5 hash for PDF/PNG figures |
| `METHODSX_PIPELINE_output_manifest.csv` | Size, timestamp, and MD5 hash for current analytical and rendered outputs |

## Interpretation boundary

The empirical application describes aggregate application preferences. It does
not establish realized student mobility, causal effects of policy changes,
individual continuity across multiple edges, or institutional quality.


# Contribution and independent validation protocol

## Scope

Changes to analysis code, raw aggregate inputs, simulation design, manuscript
claims, or figure code must pass the same staged validation before release.

## Contributor workflow

1. Work in a dedicated branch.
2. Record the scientific reason for the change.
3. Run the diagnostic gate from `REPRODUCIBILITY.md`.
4. Inspect every failed or changed validation row.
5. Run a clean publication reproduction for changes affecting analytical
   results or release outputs.
6. Compare configuration, input, figure, and output manifests.
7. Update the manuscript only from the regenerated tables.
8. Request independent review before merging.

## Second-author software validation

Kornél Dénes's planned validation should include:

- checking that the documented environment variables control the intended
  configuration;
- tracing raw aggregate inputs through annual graph construction;
- reviewing the alternating cut lower-bound implementation, both signs,
  deterministic starts, random restarts, and retained random floor;
- verifying the simulation calibration and null false-positive calculation;
- checking estimator and resolution robustness summaries;
- reproducing a diagnostic run independently;
- reproducing or auditing the publication run and its 16 acceptance checks;
- checking that numerical manuscript statements match the CSV outputs;
- visually inspecting graphical abstract, HTML, and DOCX products;
- confirming the final CRediT statement and approving the manuscript.

Validation findings should be recorded in an issue or signed release checklist.
Authorship is not considered submission-ready until the contribution has been
performed and both authors accept accountability for the final work.

## Reporting issues

An issue report should include the operating system, R version, effective
pipeline configuration, the relevant validation rows, and the smallest code or
data example that reproduces the problem. Do not attach restricted or
individual-level administrative data.


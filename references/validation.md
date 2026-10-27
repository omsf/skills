# Validation record

Checked 2026-10-26. This records bounded checks, not certification of scientific
accuracy or production readiness.

## Collection checks

- All **16** directory workflow titles in `workflows.txt` map one-to-one to existing
  `*/SKILL.md` files in the root README.
- All 16 skills pass `uvx --from skills-ref agentskills validate <skill-directory>`
  (`agentskills` version 0.1.1).
- All skills have matching directory/frontmatter names, descriptions, source
  revision permalinks, and fewer than 500 lines.
- All Python examples parse with `ast.parse`; all shell examples pass `bash -n`;
  all frontmatter and YAML examples parse with `yaml.safe_load`.
- Local Markdown links resolve. Pinned GitHub source-file links were checked.
  Some documentation sites reject automated HEAD requests; GET checks confirmed
  the OpenFF FAQ, OpenFE showcase, and corrected QCArchive dataset-guide URL.

## Targeted runtime checks

Small CPU checks used an isolated temporary environment outside the repository:
OpenFF Toolkit 0.19.0, Interchange 0.5.5, OpenMM 8.6.1, RDKit 2026.9.1,
OpenFF Units 0.4.0, NAGL 0.6.1, NAGL Models 2026.9.0.

- Executed the Python examples for conformer energies, force-field modification,
  parameter inspection, and QCArchive conversion; all completed successfully.
  QCArchive used the live public instance, qcportal 0.71, QCEngine 0.52.0, and
  qcelemental 0.51.4. Public entries may change over time.
- Exported the 27-atom example to six nonempty AMBER/GROMACS files and evaluated
  its OpenMM point energy. The expected AMBER switching-function warning appeared;
  this is not evidence of strict cross-engine equivalence.
- Initialized and minimized the virtual-site example and ran **10 CPU steps** with
  finite energy. This was not the full 10,000-step tutorial or TIP5P comparison.
- Tested the adapted clash-removal helper: removes clashing water, refuses a
  protein/cofactor-like non-solvent clash, preserves input topology and box vectors.
- Regression-checked OpenMM CSV parsing with its leading-`#` header intact.
- OpenADMET catalog research separately tested Intake catalog loading, local and
  remote Parquet reads, optional dependency requirements, and schema differences.

## Not established

No full OpenFE free-energy campaigns, ADMET model training, long simulations,
protein/PTM end-to-end runs, GPU validation, or AMBER/GROMACS energy comparisons
were performed. Source-checked instructions for these paths are not run-verified.
No convergence, predictive accuracy, or clinical/regulatory suitability is claimed.

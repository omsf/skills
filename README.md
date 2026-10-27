# OMSF workflow skills

Portable [Agent Skills](https://agentskills.io/specification) for every workflow in the
[OMSF workflow directory](https://directory.omsf.io/workflows/llms.txt): 16 skills
covering Open Force Field, Open Free Energy, and OpenADMET.

Each `<name>/SKILL.md` teaches an agent when to use the workflow, what to ask
for, how to prepare the environment, which steps to follow, and how to check the
result. These are agent instructions, not vendored software or automated simulation
pipelines. The layout follows the catalog-and-skill-directory pattern used by
[NVIDIA's skills repository](https://github.com/NVIDIA/skills); this project is not
an NVIDIA product or an upstream certification program.

## Use a skill

1. Pick a workflow below and review its instructions, especially environment and
   scientific limitations.
2. Copy that entire skill directory into your agent's supported skill location.
   For clients using project-local `.agents/skills/`, run from your target project:

   ```bash
   mkdir -p .agents/skills
   # Replace /path/to/omsf-skills with this checkout. Use a fresh destination.
   cp -R /path/to/omsf-skills/openfe-rbfe-cli .agents/skills/
   ```

   Claude Code uses `.claude/skills/`; other clients may use different locations.
   Reload or restart the client as needed. Installing a skill does **not** install
   its scientific dependencies.
3. Ask for the workflow and supply your inputs. For example:

   > Use openfe-rbfe-cli to plan relative binding free energies for my aligned
   > ligand SDF and prepared protein. Inspect the network first; do not launch
   > simulations until I approve the settings and compute budget.

Skills are independently installable. No global tool permissions or agent-specific
plugins are required by this repository.

## Workflow catalog

Workflow titles match the [OMSF workflow directory](https://directory.omsf.io/workflows/llms.txt)
as retrieved **2026-10-26**.
Each skill links to the upstream material it was researched against. Revision
pins identify source material, not a guarantee that all dependency versions remain
available or compatible.

| Directory workflow | Skill |
| --- | --- |
| Export an OpenFF OpenMM System to AMBER and GROMACS files | [openff-export-amber-gromacs](openff-export-amber-gromacs/SKILL.md) |
| Generating Partial Charges with the OpenFE CLI | [openfe-partial-charges](openfe-partial-charges/SKILL.md) |
| Integrating a Sage ligand with an Amber protein | [openff-sage-amber-protein](openff-sage-amber-protein/SKILL.md) |
| OpenADMET Data Catalogs | [openadmet-data-catalogs](openadmet-data-catalogs/SKILL.md) |
| OpenADMET Showcase | [openadmet-showcase](openadmet-showcase/SKILL.md) |
| OpenFE Showcase | [openfe-showcase](openfe-showcase/SKILL.md) |
| OpenFF Example: Conformer Energies | [openff-conformer-energies](openff-conformer-energies/SKILL.md) |
| OpenFF Example: Forcefield Modification | [openff-forcefield-modification](openff-forcefield-modification/SKILL.md) |
| OpenFF Example: Inspecting SMIRNOFF Parameters | [openff-inspect-parameters](openff-inspect-parameters/SKILL.md) |
| OpenFF Example: QCArchive Interface | [openff-qcarchive](openff-qcarchive/SKILL.md) |
| OpenFF Example: SMIRNOFF Simulation | [openff-smirnoff-simulation](openff-smirnoff-simulation/SKILL.md) |
| OpenFF Example: Virtual Sites Showcase | [openff-virtual-sites](openff-virtual-sites/SKILL.md) |
| OpenFF Toolkit Showcase | [openff-toolkit-showcase](openff-toolkit-showcase/SKILL.md) |
| Post Translational Modification Prototype | [openff-ptm-prototype](openff-ptm-prototype/SKILL.md) |
| Relative Free Energies with the OpenFE CLI | [openfe-rbfe-cli](openfe-rbfe-cli/SKILL.md) |
| Visualizations with OpenFF | [openff-visualization](openff-visualization/SKILL.md) |

## Limitations

A successful example run is not evidence of convergence, predictive accuracy,
or suitability for clinical or regulatory decisions. Full end-to-end scientific
validation is not claimed.

Repository maintenance and agent execution instructions are in
[AGENTS.md](AGENTS.md).

## Attribution and licensing

Sources and upstream license identifiers are recorded in individual skills.
Those identifiers describe upstream projects, not a grant
of rights to this collection or every underlying dataset. Check dataset-specific
terms separately. No repository-wide license has been selected for this new
collection; select one before public distribution. Upstream notebooks and datasets
are linked rather than bundled.

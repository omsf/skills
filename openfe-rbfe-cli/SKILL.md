---
name: openfe-rbfe-cli
description: Plan, execute, and gather relative binding or hydration free energy calculations with the OpenFE CLI (`openfe plan-rbfe-network`, `openfe plan-rhfe-network`, `openfe quickrun`, `openfe gather`) without writing Python. Use when a user wants an alchemical free-energy campaign from SDF ligands plus a protein structure, needs to know the network_setup/transformations JSON layout, needs to submit quickrun legs to SLURM/HPC, or needs to turn result JSONs into a DG/DDG table.
---

# Relative Free Energies with the OpenFE CLI

## Purpose and scope

Three CLI commands cover an RBFE/RHFE campaign:

1. `openfe plan-rbfe-network` — define systems, generate the ligand network and
   one transformation JSON per alchemical leg.
2. `openfe quickrun` — execute one leg, with the repeat count stored in its JSON.
3. `openfe gather` — combine result JSONs into a TSV of `DG` or `DDG`.

`plan-rhfe-network` is the hydration analogue: same run/gather steps, no
protein. Scope: RBFE/RHFE via the CLI. ABFE, SepTop, and Python-level
customisation are out of scope.

## Required inputs and clarifications

Ask before starting:

- **Ligands**: SDF or MOL2 file, or a directory of them. Molecule names from the
  SDF become ligand identifiers and appear in every output filename.
- **Protein**: a PDB or PDBx/mmCIF file. For membrane systems use
  `--protein-membrane` with a *fully built, solvated, pre-equilibrated* system
  including correct box vectors — not a bare protein PDB.
- **Repeats and budget.** Default is 3 independent sampling repeats per leg and
  5 ns production per repeat. Confirm the compute budget.
- **Charge set.** If you pre-generated charges with `openfe charge-molecules`
  (see `openfe-partial-charges`), pass the *charged* SDF here so planner and
  legs share one charge set.
- **Approval for GPU/large runs.** Each leg is hours of GPU time. Get explicit
  user approval before launching a full network, and never run the whole network
  on a login node.

## Environment

```bash
micromamba create -c conda-forge -n openfe openfe=1.13
micromamba activate openfe
openfe --version
openfe plan-rbfe-network --help
```

The install targets OpenFE release **v1.13.0** (2026-10-08). Commands were
source-checked on `main` `16a5feaa0a7e19d2348d9d63ee15ea51fd471e15`, not executed.
Main-only options (notably newer orchestration/membrane/mmCIF features) may differ
from the release; inspect installed `--help` before use. The CLI's `plan-rbfe-network` /
`quickrun` / `gather` trio has been stable since 0.7; the exact options below
were read from that revision.

## Procedure

### 1. Fetch tutorial data (optional)

```bash
mkdir -p rbfe_work && cd rbfe_work
openfe fetch rbfe-tutorial          # tyk2_ligands.sdf, tyk2_protein.pdb, cli_tutorial.md, rbfe_python_tutorial.ipynb
openfe fetch rbfe-tutorial-results  # precomputed results
tar xzf rbfe_results.tar.gz         # -> results/
```

`openfe fetch <name> -d <dir>` chooses the output directory; the default is the
current directory. Both fetchers need internet access.

### 2. Plan the network

```bash
openfe plan-rbfe-network \
    -M tyk2_ligands.sdf \
    -p tyk2_protein.pdb \
    -o network_setup \
    --n-protocol-repeats 1 \
    -n 4
```

Options (verified flags):

| Flag | Meaning |
| --- | --- |
| `-M`, `--molecules` | ligand SDF/MOL2 file or directory (required) |
| `-p`, `--protein` | protein PDB/mmCIF; mutually exclusive with `--protein-membrane` |
| `--protein-membrane` | pre-solvated membrane PDB/mmCIF |
| `-C`, `--cofactors` | cofactor SDF, repeatable |
| `-o`, `--output-dir` | output dir, default `./alchemicalNetwork` |
| `--n-protocol-repeats` | repeats baked into each transformation JSON, default 3 |
| `-n`, `--n-cores` | processes for partial-charge generation |
| `-s`, `--settings` | YAML overriding mapper/network/partial charge |
| `--overwrite-charges` | recompute charges even if present |
| `--networks-only` | write only the AlchemicalNetwork JSON + ligand graphml |

Exactly one of `--protein` / `--protein-membrane` must be given, otherwise the
command exits with a `UsageError`.

Defaults chosen for you (from the CLI docstring at the verified revision):
Kartograf atom mapper (`atom_max_distance=0.95`, `atom_map_hydrogens=True`,
`map_hydrogens_on_hydrogens_only=True`, `map_exact_ring_matches_only=True`,
`allow_partial_fused_rings=True`, `allow_bond_breaks=False`), minimal spanning
tree network with the LOMAP default score, water with NaCl at 0.15 M
(neutralised), and the OpenMM hybrid-topology RBFE protocol.

Planning is slow because it (a) runs a semi-empirical QM job per ligand to make
AM1-BCC charges and (b) enumerates and scores all ligand pairs before extracting
the MST — hence `-n`.

### 3. (Optional) customise via YAML

```yaml
# settings.yaml
mapper:
  method: lomap                       # lomap | kartograf (class names also accepted)
network:
  method: generate_maximal_network    # see list below
partial_charge:
  method: nagl
  settings:
    nagl_model: null                  # null => latest am1bcc NAGL model
```

```bash
openfe plan-rbfe-network -M tyk2_ligands.sdf -p tyk2_protein.pdb \
    -o network_setup --n-protocol-repeats 1 -s settings.yaml
```

Accepted `network.method` values: `generate_minimal_spanning_network` (default),
`generate_minimal_redundant_network`, `generate_radial_network`,
`generate_maximal_network`, `generate_lomap_network` (aliases `mst`, `radial`).
`mapper.method` strings and `network.method` strings are lower-cased and matched
case-insensitively. Unknown top-level YAML sections are ignored with a
`RuntimeWarning` rather than an error, so a typo'd key silently does nothing —
diff your YAML against `openfe plan-rbfe-network --help`.

Projected output layout:

```text
network_setup
├── ligand_network.graphml
├── network_setup.json
└── transformations
    ├── rbfe_lig_ejm_31_complex_lig_ejm_42_complex.json
    ├── rbfe_lig_ejm_31_solvent_lig_ejm_42_solvent.json
    └── ...
```

Each RBFE comparison needs two legs (solvent + complex); `openfe gather` later
subtracts them. One JSON = one leg, not one ligand pair.

### 4. Inspect the network

```bash
openfe view-ligand-network network_setup/ligand_network.graphml
```

Caveat: the current implementation loads the graphml, calls
`plot_atommapping_network`, forces the `TkAgg` matplotlib backend, and shows a
static figure. There is no click-to-inspect-edge interactivity and it needs a
display, so it fails on headless nodes. The tutorial's "click on the edge to see
the mapping" description no longer matches the code. For headless use, load the
graphml in Python or read `ligand_network.graphml` directly.

### 5. Run the legs

One leg, one repeat, one command:

```bash
openfe quickrun \
    network_setup/transformations/rbfe_lig_ejm_31_solvent_lig_ejm_42_solvent.json \
    -o results_0/rbfe_lig_ejm_31_solvent_lig_ejm_42_solvent.json \
    -d results_0/rbfe_lig_ejm_31_solvent_lig_ejm_42_solvent
```

- `-o` (output option) is the result JSON path; `-d` / `--work-dir` is the
  simulation working directory (created if missing). Default work dir is the
  current directory.
- `--resume` re-uses cached task outputs after an interruption; the command
  prints the exact resume command.
- Use a **unique** output JSON and work dir per leg *and* per repeat, otherwise
  results overwrite each other.

If you ran with `--n-protocol-repeats 1`, execute each leg at least twice
(three times is the recommendation), into different result directories.

Generate job scripts only; do not submit them yet. Derive paths from the basename
(the tutorial's `${file:30}` slicing breaks when the directory name changes). Add
your site's GPU, memory, walltime, partition/account and environment directives,
review the scripts and obtain approval before invoking `sbatch`:

```bash
mkdir -p slurm
for file in network_setup/transformations/*.json; do
  leg=$(basename "$file" .json)
  for repeat in 0 1 2; do
    jobpath="slurm/${leg}_${repeat}.job"
    {
      echo '#!/usr/bin/env bash'
      echo '#SBATCH --job-name=openfe-rbfe'
      echo 'set -euo pipefail'
      printf '%q ' micromamba run -n openfe openfe quickrun "$file" \
        -o "results_${repeat}/${leg}.json" -d "results_${repeat}/${leg}"
      printf '\n'
    } > "$jobpath"
    # After adding site-specific directives and getting approval: sbatch "$jobpath"
  done
done
```

### 6. Gather results

```bash
openfe gather results/ --report dg -o final_results.tsv
```

Multiple directories are accepted: `openfe gather results_0/ results_1/ ...`.
Directories are walked recursively and any file ending in `.json` is collected;
malformed files are ignored.

| `--report` | Columns |
| --- | --- |
| `dg` (default) | `ligand`, `DG(MLE) (kcal/mol)`, `uncertainty (kcal/mol)` |
| `ddg` | `ligand_i`, `ligand_j`, `DDG(i->j) (kcal/mol)`, `uncertainty (kcal/mol)` |
| `raw` | `leg`, `ligand_i`, `ligand_j`, `DG(i->j) (kcal/mol)`, `MBAR uncertainty (kcal/mol)` |

Other useful flags: `-o -` / default writes to stdout; `--tsv` keeps stdout
tab-separated; `--allow-partial` warns-and-skips instead of exiting on
incomplete edges.

## Validation and expected outputs

```bash
# A sane dg table: 10 ligands, ~0-2 kcal/mol spread, uncertainties ~0.05-0.2
column -t final_results.tsv | head
```

Stored tutorial output (TYK2, 10 ligands; not a fresh simulation). `dg` values
are on an arbitrary reference scale, not physical absolute binding free energies:

```text
ligand	DG(MLE) (kcal/mol)	uncertainty (kcal/mol)
lig_ejm_31	-0.09	0.05
lig_ejm_42	0.7	0.1
lig_ejm_46	-0.98	0.05
...
```

Sanity checks before trusting numbers: `--report raw` should show both `solvent`
and `complex` legs for every edge; large uncertainties warrant checking sampling and mapping. No single uncertainty
threshold establishes convergence.

## Failure modes and safety

- **`UsageError: Either --protein or --protein-membrane must be provided.`** —
  supply one; not both.
- **`ERROR: Every edge must have at least two simulation repeats`** — the `dg`
  report needs ≥2 repeats per leg. With `--n-protocol-repeats 1` you must run
  `quickrun` at least twice per leg. `--report ddg` and `--report raw` do not
  have this constraint.
- **`ERROR: The results network is disconnected due to failed or missing edges`**
  — a common relative-affinity scale cannot be inferred for disconnected components. Fix or add edges,
  or use `--report ddg`, or force `--allow-partial` (which may bias the MLE).
- **Silent per-leg failures** — a failing leg does not abort the campaign; the
  failure is recorded inside the result JSON. Grep the JSONs for failed units
  rather than assuming success from file existence.
- **Overwritten results** — `quickrun` happily writes over existing outputs.
  Never reuse one `-o`/`-d` pair across repeats.
- **`view-ligand-network` needs a display** (TkAgg). Headless = crash.
- **Cost safety**: an RBFE leg is GPU-hours. Confirm with the user before
  planning/running a full network, before raising `--n-protocol-repeats`, and
  before changing lambda windows or production length upward.
- **Newer orchestration exists.** `openfe setup-task-campaign`, `openfe
  run-task`, `openfe status`, and `openfe update-task-db` (OpenFE ≥1.12/1.13)
  wrap the JSON/warehouse loop into a task database for large campaigns. The
  tutorial predates them; consider them for campaigns bigger than a few edges,
  but verify they exist in your installed version first.
- **Do not "fix" settings silently.** `--n-protocol-repeats`, charge method, and
  lambda/production lengths are scientific choices. Report defaults, do not
  quietly upgrade them.

## Sources

Tutorial (rendered docs, checked 2026-10-26):
<https://docs.openfree.energy/en/latest/tutorials/rbfe_cli_tutorial.html>

Upstream tutorial source, OpenFreeEnergy/ExampleNotebooks @
`7c1818c6d8954daf3df902fa6e5e224226dc0a5d` (main, 2026-10-26):
- <https://github.com/OpenFreeEnergy/ExampleNotebooks/blob/7c1818c6d8954daf3df902fa6e5e224226dc0a5d/rbfe_tutorial/cli_tutorial.md>
- <https://github.com/OpenFreeEnergy/ExampleNotebooks/blob/7c1818c6d8954daf3df902fa6e5e224226dc0a5d/rbfe_tutorial/README.md>
- Directory: <https://github.com/OpenFreeEnergy/ExampleNotebooks/tree/main/rbfe_tutorial>
  (tutorial data `tyk2_ligands.sdf` / `tyk2_protein.pdb` are byte-identical to
  the copies under `showcase/inputs/`)

CLI implementation read, OpenFreeEnergy/openfe @
`16a5feaa0a7e19d2348d9d63ee15ea51fd471e15` (main, 2026-10-26); release
`v1.13.0`:
- <https://github.com/OpenFreeEnergy/openfe/blob/16a5feaa0a7e19d2348d9d63ee15ea51fd471e15/src/openfecli/commands/plan_rbfe_network.py>
- <https://github.com/OpenFreeEnergy/openfe/blob/16a5feaa0a7e19d2348d9d63ee15ea51fd471e15/src/openfecli/commands/quickrun.py>
- <https://github.com/OpenFreeEnergy/openfe/blob/16a5feaa0a7e19d2348d9d63ee15ea51fd471e15/src/openfecli/commands/gather.py>
- <https://github.com/OpenFreeEnergy/openfe/blob/16a5feaa0a7e19d2348d9d63ee15ea51fd471e15/src/openfecli/commands/view_ligand_network.py>
- <https://github.com/OpenFreeEnergy/openfe/blob/16a5feaa0a7e19d2348d9d63ee15ea51fd471e15/src/openfecli/parameters/plan_network_options.py>
- <https://github.com/OpenFreeEnergy/openfe/blob/16a5feaa0a7e19d2348d9d63ee15ea51fd471e15/src/openfecli/fetchables.py>
- Protocol defaults: <https://github.com/OpenFreeEnergy/openfe/blob/16a5feaa0a7e19d2348d9d63ee15ea51fd471e15/src/openfe/protocols/openmm_rfe/hybridtop_protocols.py>
- Install docs: <https://docs.openfree.energy/en/latest/installation.html>

Licensing: OpenFE and the ExampleNotebooks are MIT; this skill is an original
summary, not a copy of upstream text.
---
name: openfe-partial-charges
description: Generate and store partial charges for a set of small molecules ahead of OpenFE free-energy planning, using the `openfe charge-molecules` CLI. Use when preparing ligands (or cofactors) for OpenFE RBFE/RHFE/AHFE campaigns, when you need one reproducible charge set shared across transformation legs, when picking a charge method (am1bcc, am1bccelf10, nagl, espaloma, forcefield), or when a plan-rbfe-network run spends a long time generating AM1-BCC charges repeatedly.
---

# Generating Partial Charges with the OpenFE CLI

## Purpose and scope

`openfe charge-molecules` runs partial-charge generation once, for a whole ligand
set, and writes the charges into the `atom.dprop.PartialCharge` SD tag of an SDF.
Downstream planners (`openfe plan-rbfe-network`, `plan-rhfe-network`) then reuse
those charges instead of recomputing them for every leg.

Why this matters: the OpenFE planners assign charges while building the network,
so re-planning the same ligands several times to search for a better network
re-pays the cost each time. More importantly, recomputing charges per leg can
give the same ligand different charges in the solvent and complex legs, which
breaks the thermodynamic cycle (see the tutorial's Osato et al. reference).

Scope: charge generation only. This skill does not plan or run simulations.
For campaign planning/execution see `openfe-rbfe-cli`.

## Required inputs and clarifications

Ask before starting if any of these are unclear:

- **Input file(s)**: one `.sdf`, one `.mol2`, or a directory. Directories are
  scanned non-recursively for `.sdf` and `.mol2` only.
- **Which charge method?** Default `am1bcc` on the AmberTools backend is the
  documented OpenFE CLI default; the selected force field may instead specify NAGL.
  Confirm the charge method together with the force field. Do not silently swap methods:
  `am1bccelf10`, `nagl`, `espaloma`, and `forcefield` are different scientific
  choices with different accuracy/cost trade-offs.
- **Extra packages**: `am1bccelf10` needs OpenEye; `nagl` needs `openff-nagl`;
  `espaloma` needs `espaloma_charge`; `forcefield` needs explicit offxml files.
  Confirm these are installed before promising a method.
- **Approval**: local charge backends do not inherently upload molecules. Ask before
  any remote service, installing optional packages, or downloading a NAGL model. For large sets (> a few hundred molecules) on a
  shared machine, confirm CPU/time budget with the user.

## Environment

OpenFE ships both the core library and the CLI, and is POSIX-only (macOS/Linux).
Recommended install is micromamba from conda-forge. The commands below are
source-checked against the main revision cited in Sources, not executed here; check
your installed `--help`, because main-only options may not be in a release yet:

```bash
micromamba create -c conda-forge -n openfe openfe=1.13
micromamba activate openfe
openfe --version
openfe charge-molecules --help
```

Alternative lock-file installation (needs internet): the `latest` URL moves, so
record the resolved release and lock-file checksum before claiming reproducibility:

```bash
curl -LOJ https://github.com/OpenFreeEnergy/openfe/releases/latest/download/openfe-conda-lock.yml
micromamba create -n openfe --file openfe-conda-lock.yml
micromamba activate openfe
```

Version caveats verified against this source revision:

- `charge-molecules` was added in OpenFE 1.3 (`requires_ofe=(1, 3)`); it does
  **not** exist in older 0.x installs. Latest release checked: **v1.13.0**
  (2026-10-08).
- The command and its help text live in
  `openfecli/commands/generate_partial_charges.py`; the YAML keys are validated
  by `openfecli/parameters/plan_network_options.py`.
- `am1bcc` requires the AmberTools backend (default). The `rdkit` backend is
  rejected for `am1bcc`/`am1bccelf10`.

## Procedure

### 1. Get input molecules

If you have no data yet, the OpenFE RBFE tutorial ligands are fetchable:

```bash
mkdir -p charge_work && cd charge_work
openfe fetch rbfe-tutorial
ls   # tyk2_ligands.sdf  tyk2_protein.pdb  cli_tutorial.md  rbfe_python_tutorial.ipynb
```

Inspect molecule names and count before charging, because names become the
ligand identifiers used by every downstream filename:

```bash
python -c "
from rdkit import Chem
mols = [m for m in Chem.SDMolSupplier('tyk2_ligands.sdf', removeHs=False)]
print(len(mols), [m.GetProp('_Name') for m in mols])
"
```

### 2. Generate charges

```bash
openfe charge-molecules -M tyk2_ligands.sdf -o charged_tyk2_ligands.sdf
```

- `-M` / `--molecules` accepts a single SDF/MOL2 file or a directory.
- `-o` / `--output` is the output SDF path. **The output file must not already
  exist**; the command raises `FileExistsError` rather than overwriting. Pick a
  new name or delete the old file first.

Parallelise with `-n` / `--n-cores` (AM1-BCC is a semi-empirical QM step, so
this is where the wall-clock time goes). Ordering is preserved in the output:

```bash
openfe charge-molecules -M tyk2_ligands.sdf -o charged_tyk2_ligands.sdf -n 4
```

### 3. (Optional) customise the charge method

Write a settings YAML. Only the `partial_charge` section is honoured by
`charge-molecules`; `mapper` and `network` sections are ignored.

```yaml
# settings.yaml
partial_charge:
  method: am1bcc          # am1bcc | am1bccelf10 | nagl | espaloma | forcefield
  settings:
    off_toolkit_backend: ambertools   # ambertools | openeye | rdkit
    number_of_conformers: null        # null => reuse the input conformer
```

```bash
openfe charge-molecules -M tyk2_ligands.sdf -o charged_tyk2_ligands.sdf \
    -n 4 -s settings.yaml
```

Full documented option set (some need extra packages):

```yaml
partial_charge:
  method: am1bcc
  # method: am1bccelf10     # requires off_toolkit_backend: openeye; use ~500 conformers
  # method: nagl            # requires openff-nagl
  # method: espaloma        # requires espaloma_charge; ambertools or rdkit backend only
  # method: forcefield      # requires `forcefields` below
  settings:
    off_toolkit_backend: ambertools
    number_of_conformers: null
    # nagl_model: null                    # null => latest am1bcc NAGL model
    # forcefields: [openff-2.3.0.offxml]  # used by `method: forcefield`
```

Method names are matched case-insensitively (values are lower-cased), so
`method: NAGL` and `method: nagl` are equivalent.

### 4. Reuse the charges downstream

Pass the charged SDF straight to the planner instead of the raw one:

```bash
openfe plan-rbfe-network -M charged_tyk2_ligands.sdf -p tyk2_protein.pdb \
    -o network_setup --n-protocol-repeats 3 -n 4
```

By default `charge-molecules` **only charges molecules that have no charges**;
ligands that already carry a `atom.dprop.PartialCharge` tag are passed through
untouched. Zeroed charges are always recomputed. To force new charges, add
`--overwrite-charges`:

```bash
openfe charge-molecules -M charged_tyk2_ligands.sdf -o recharged.sdf \
    --overwrite-charges -s settings.yaml
```

The same `--overwrite-charges` flag exists on `plan-rbfe-network` /
`plan-rhfe-network`; without it they likewise respect pre-existing charges.

## Validation and expected outputs

The CLI prints a banner and a summary; check that the method and molecule list
match what you asked for:

```text
SMALL MOLECULE PARTIAL CHARGE GENERATOR
_________________________________________

Parsing in Files:
	Got input:
		Small Molecules: SmallMoleculeComponent(name=lig_ejm_31) ...
Using Options:
	Partial Charge Generation: am1bcc
```

Then verify the output file actually carries charges:

```bash
grep -c 'atom.dprop.PartialCharge' charged_tyk2_ligands.sdf   # == number of ligands
```

```python
# Check the charge vector length matches atom count, per molecule
from rdkit import Chem
import math
count = 0
for mol in Chem.SDMolSupplier("charged_tyk2_ligands.sdf", removeHs=False):
    assert mol is not None, "invalid SDF record"
    v = [float(x) for x in mol.GetProp("atom.dprop.PartialCharge").split()]
    assert all(math.isfinite(q) for q in v), "non-finite charge"
    assert abs(sum(v) - Chem.GetFormalCharge(mol)) < 1e-3, "net charge mismatch"
    assert len(v) == mol.GetNumAtoms(), (mol.GetProp("_Name"), len(v), mol.GetNumAtoms())
    count += 1
assert count > 0, "no molecules read"
print("ok:", count, "molecules charged")
```

Reproducibility check: re-charge the same input with the same settings into a
fresh file and compare charge vectors, not file bytes (SDF writers can reorder
whitespace).

## Failure modes and safety

- **`FileExistsError: The output file ... already exists`** — intentional.
  Choose a new output name or delete the existing file.
- **`Unable to find molecules in <path>`** — wrong extension or an empty
  directory. Directory loading is not recursive and ignores other formats.
- **`ExceptionGroup` / per-molecule errors** — `charge-molecules` calls
  `bulk_assign_partial_charges(..., raise_errors=True)` on purpose, so one bad
  molecule aborts the run with its name and SMILES attached as a note. Fix the
  input structure (e.g. missing hydrogens/explicit conformer) and rerun.
- **`generate_n_conformers` errors** — a value >1 is only valid for
  `am1bccelf10`. `am1bcc`, `nagl`, and `espaloma` allow at most one conformer.
  Also, if `number_of_conformers` is `null` the input molecule must already
  have a conformer.
- **Missing backend packages** — OpenEye / NAGL / espaloma charge methods fail
  at charge time, not at parse time. Verify the package is importable in the
  same environment before launching a long run.
- **MOL2 input** — the CLI loads at most one molecule per MOL2 file; prefer SDF
  for multi-molecule sets.
- **Do not change the charge method mid-campaign.** Mixing charge sets between
  legs or between campaign revisions silently invalidates the cycle closure
  assumption. If you must change methods, re-charge and redo the campaign.
- Resource safety: `-n` spawns that many worker processes, each running AM1-BCC.
  Do not set `-n` above the cores you actually own; ask the user before large
  batches.

## Sources

Tutorial (rendered docs, checked 2026-10-26):
<https://docs.openfree.energy/en/latest/tutorials/charge_molecules_cli_tutorial.html>

Upstream tutorial source, OpenFreeEnergy/ExampleNotebooks @
`7c1818c6d8954daf3df902fa6e5e224226dc0a5d` (main, 2026-10-26):
- <https://github.com/OpenFreeEnergy/ExampleNotebooks/blob/7c1818c6d8954daf3df902fa6e5e224226dc0a5d/cli_tutorials/cli_charge_molecules.md>
- Directory: <https://github.com/OpenFreeEnergy/ExampleNotebooks/tree/main/cli_tutorials>
- Environment pin reference: <https://github.com/OpenFreeEnergy/ExampleNotebooks/blob/7c1818c6d8954daf3df902fa6e5e224226dc0a5d/environment.yaml>
  (note: that file pins `openfe==1.13.*` from conda-forge *and* pip-installs
  `git+https://github.com/OpenFreeEnergy/openfe@main`; the latter wins, so
  notebook CI tracks `main`, not the 1.13 pin)

CLI implementation read, OpenFreeEnergy/openfe @
`16a5feaa0a7e19d2348d9d63ee15ea51fd471e15` (main, 2026-10-26); release
`v1.13.0` (2026-10-08):
- <https://github.com/OpenFreeEnergy/openfe/blob/16a5feaa0a7e19d2348d9d63ee15ea51fd471e15/src/openfecli/commands/generate_partial_charges.py>
- <https://github.com/OpenFreeEnergy/openfe/blob/16a5feaa0a7e19d2348d9d63ee15ea51fd471e15/src/openfecli/parameters/misc.py>
- <https://github.com/OpenFreeEnergy/openfe/blob/16a5feaa0a7e19d2348d9d63ee15ea51fd471e15/src/openfecli/parameters/molecules.py>
- <https://github.com/OpenFreeEnergy/openfe/blob/16a5feaa0a7e19d2348d9d63ee15ea51fd471e15/src/openfe/protocols/openmm_utils/charge_generation.py>
- <https://github.com/OpenFreeEnergy/openfe/blob/16a5feaa0a7e19d2348d9d63ee15ea51fd471e15/src/openfe/protocols/openmm_utils/omm_settings.py>
- Install docs: <https://docs.openfree.energy/en/latest/installation.html>

Licensing: OpenFE and the ExampleNotebooks are MIT; this skill is an original
summary, not a copy of upstream text.

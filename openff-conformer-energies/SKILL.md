---
name: openff-conformer-energies
description: Compute vacuum single-point and minimized conformer energies for a small molecule with a SMIRNOFF/Sage force field via the OpenFF Toolkit and Interchange, for the openff-toolkit conformer_energies example. Use for requests about ranking SDF conformers by energy, energy minimization of small molecules, RMSD between input and minimized conformers, or the conformer_energies.py CLI.
---

# OpenFF Example: Conformer Energies

## Purpose and scope

Read several conformers of **one** chemical species from an SDF file, apply a
SMIRNOFF force field, and report for each conformer the vacuum potential energy
before and after minimization plus the RMSD between the two geometries. Writes one
minimized SDF per conformer and a CSV summary.

Scope ends at vacuum energies of a single small molecule. It is *not* a solvated
free-energy method, and absolute vacuum energies are only comparable within one
force field / charge model (relative conformer energies are the meaningful output).

Upstream files: `examples/conformer_energies/conformer_energies.ipynb`,
`examples/conformer_energies/conformer_energies.py` (CLI),
`examples/conformer_energies/README.md`,
`examples/conformer_energies/ruxolitinib_conformers.sdf`.

## Required inputs and clarifications

1. **The molecule**: an SDF whose records are all the *same* species (conformers).
   The notebook ships `ruxolitinib_conformers.sdf` (10 records) and reaches it via
   `get_data_file_path("molecules/ruxolitinib_conformers.sdf")`.
2. **Force field**: the current example hardcodes `openff_unconstrained-2.3.0.offxml`
   (Sage 2.3.0). **Do not silently switch versions or switch to the constrained file.**
   The unconstrained file is deliberate — minimization of hydrogen bonds with
   constraints in place gives wrong energies (see the toolkit FAQ on "unconstrained").
   Note the upstream `README.md` still says the CLI uses "Sage (openff-2.1.0)"; the
   code says 2.3.0. Trust the code, and confirm with the user if they need 2.1.0.
3. **Charge model**: Sage 2.3.0 specifies `<NAGLCharges model_file="openff-gnn-am1bcc-1.0.0.pt">`,
   i.e. ML (NAGL) AM1-BCC-like charges, **not** an AmberTools/OpenEye AM1-BCC call.
   So the notebook's warning about AM1-BCC implementation differences only applies if
   you deliberately downgrade to a force field that uses `ToolkitAM1BCC`.
4. **How many conformers**: 10 ruxolitinib conformers takes ~15 s on CPU. For hundreds
   of conformers or larger molecules, ask before running; charge assignment dominates.
5. **GPU**: none is needed or used by this example. Do not add GPU work unasked.

## Environment setup and version caveats

```bash
mamba create -n openff-toolkit -c conda-forge "openff-toolkit=0.19.0"
mamba activate openff-toolkit
python -c "from openff.toolkit import GLOBAL_TOOLKIT_REGISTRY; print(GLOBAL_TOOLKIT_REGISTRY.registered_toolkit_versions)"
# e.g. {'OpenFF NAGL': '0.6.1', 'The RDKit': '2026.9.1', 'Built-in Toolkit': None}
```

- The canonical environment for running the notebooks is
  `devtools/conda-envs/examples.yaml` (name `openff-toolkit-examples`); CI runs
  `pytest --nbval-lax --dist loadscope --ignore=examples/deprecated examples`.
- **You need `openff-nagl` + `openff-nagl-models`** (pulled in by `openff-toolkit`) for
  Sage 2.3.0 charges. They drag in `libtorch`, so the solve is large.
- The Toolkit supplies `ForceField.create_interchange`; Interchange supplies
  `Interchange.positions`, `Interchange.minimize`, and
  `openff.interchange.drivers.openmm.get_openmm_energies`. `openff-toolkit` alone is
  not enough.
- `get_data_file_path` resolves inside the *installed* `openff-toolkit` package; the
  example SDF is not in your working directory.
- Verified against `openff-toolkit` 0.19.0, `openff-interchange` 0.5.5, `openff-units`
  0.4.0, `openff-nagl` 0.6.1, `openff-nagl-models` 2026.9.0, RDKit 2026.9.1,
  OpenMM 8.6.1, and `openff-2.3.0.offxml` sha256
  `7a0d7195a4b717e29fe14aa2cf1bc3203e87ce4aef9cc2ef10fb5145b7edcf63`.

## Steps

Run in a fresh scratch directory; the example writes SDF and CSV files.

### 1. Load and collapse the SDF

`Molecule.from_file` returns a `Molecule` for a one-record file and a `list` for a
multi-record file, so normalize both cases:

```python
from openff.toolkit import ForceField, Molecule
from openff.toolkit.utils import get_data_file_path

sdf_path = get_data_file_path("molecules/ruxolitinib_conformers.sdf")
loaded = Molecule.from_file(sdf_path)

# from_file returns a Molecule for a single-record file, a list otherwise
if isinstance(loaded, Molecule):
    loaded = [loaded]

# Each SDF record becomes its own Molecule, so collapse them into one
molecule = loaded.pop(0)
for next_molecule in loaded:
    if next_molecule == molecule:
        for conformer in next_molecule.conformers:
            molecule.add_conformer(conformer)
    else:
        raise ValueError("Multiple chemical species loaded")

if not molecule.name:
    molecule.name = molecule.to_hill_formula()

# Names become paths below: keep them in this working directory.
import re
molecule.name = re.sub(r"[^A-Za-z0-9._-]+", "_", molecule.name).strip("._") or "molecule"

print(f"Loaded {molecule.n_conformers} conformers of {molecule.name}")
```

### 2. Parametrize once, then reuse

Charges are the expensive part; build one `Interchange` and move its positions per
conformer instead of re-parametrizing:

```python
from openff.interchange.drivers.openmm import get_openmm_energies

forcefield = ForceField("openff_unconstrained-2.3.0.offxml")
interchange = forcefield.create_interchange(molecule.to_topology())

initial_energies, minimized_energies = [], []
minimized_molecule = Molecule(molecule)
minimized_molecule.conformers.clear()

for conformer in molecule.conformers:
    interchange.positions = conformer
    initial_energies.append(get_openmm_energies(interchange).total_energy)
    interchange.minimize(engine="openmm")
    minimized_energies.append(get_openmm_energies(interchange).total_energy)
    minimized_molecule.add_conformer(interchange.positions.to("angstrom"))
```

`get_openmm_energies(...)` returns an `EnergyReport` whose only stable per-term API in
`openff-interchange` 0.5.x is `.total_energy` plus `.energies` (a dict keyed `Bond`,
`Angle`, `Torsion`, `Nonbonded`). Older attribute names such as `.bond_energy` no
longer exist.

`Interchange.minimize(engine="openmm")` updates `interchange.positions` **in place**,
builds a throwaway `Simulation` internally (a `LangevinMiddleIntegrator` at 293.15 K is
used only as a container), and stops at a **force** tolerance of
`10 kJ/mol/nm` after at most `10_000` iterations. Both are overridable:
`interchange.minimize(engine="openmm", force_tolerance=1 * unit.kilojoule_per_mole / unit.nanometer)`.
Tighten the tolerance if you need minimized geometries to be comparable across runs.

### 3. Report energies, RMSD and files

```python
from rdkit.Chem import rdMolAlign

n_confs = molecule.n_conformers
working_mol = Molecule(molecule)
rows = [["Conformer", "Initial PE (kcal/mol)", "Minimized PE (kcal/mol)",
         "RMSD (Angstrom)"]]

for i, (init_e, init_xyz, min_e, min_xyz) in enumerate(
    zip(initial_energies, molecule.conformers, minimized_energies,
        minimized_molecule.conformers), start=1
):
    working_mol.conformers.clear()
    working_mol.add_conformer(min_xyz)
    working_mol.to_file(f"{molecule.name}_conf{i}_minimized.sdf", file_format="sdf")

    working_mol.add_conformer(init_xyz)
    rdmol = working_mol.to_rdkit()
    rmslist = []
    rdMolAlign.AlignMolConformers(rdmol, RMSlist=rmslist)

    rows.append([i, init_e.m_as("kilocalories_per_mole"),
                 min_e.m_as("kilocalories_per_mole"), rmslist[0]])
    print(f"{i:5d} / {n_confs:5d} :  {rows[-1][1]:8.3f} kcal/mol "
          f"{rows[-1][2]:8.3f} kcal/mol {rows[-1][3]:8.3f} A")

with open(f"{molecule.name}.csv", "w") as handle:
    handle.write(", ".join(str(x) for x in rows.pop(0)) + "\n")
    for line in rows:
        handle.write("{}, {:.3f}, {:.3f}, {:.3f}\n".format(*line))
```

### 4. Alternative: the shipped CLI

Download the script first (it is not installed as a command by the Toolkit):

```bash
curl --fail --location --output conformer_energies.py \
  https://raw.githubusercontent.com/openforcefield/openff-toolkit/ef074d58a2e7b03625aae7cd3aff87a8ee46c38b/examples/conformer_energies/conformer_energies.py
python conformer_energies.py -f ruxolitinib_conformers.sdf
```

**Caveat (verified):** the `-f/--filename` argument is parsed but never used — the
script always loads `get_data_file_path("molecules/ruxolitinib_conformers.sdf")`.
To use your own file, copy the script and replace that line, or use the notebook
cells above. The script also always uses `openff_unconstrained-2.3.0.offxml`.

## Validation and expected outputs

- Console table with one line per conformer, then `ruxolitinib_conf<N>_minimized.sdf`
  (one per conformer) and `ruxolitinib.csv`. The example directory's `.gitignore`
  lists exactly `ruxolitinib.csv` and `ruxolitinib_conf*_minimized.sdf`.
- Reproduced values for the shipped ruxolitinib data (kcal/mol, Sage 2.3.0, NAGL
  charges): initial energies span about `-95.5` to `-50.8`; minimized energies cluster
  tightly at about `-125.0` to `-122.4`; RMSDs 0.28-0.94 A. Every minimized energy must
  be no higher than its initial energy within numerical tolerance (an already-minimized
  input can be unchanged). The narrowing observed here is specific to this example.
- **Minimization is not bit-reproducible.** Two back-to-back runs of the identical
  script differed by up to ~0.03 A in RMSD and ~0.1 kcal/mol in minimized energy (the
  default 10 kJ/mol/nm force tolerance stops at a slightly different point near the
  minimum). Report observed precision and the minimizer tolerance; do not infer a
  universal percentage error from this one example.
- The notebook's own caveat: absolute vacuum energies are sensitive to small partial
  charge changes (e.g. different AM1-BCC implementations), while *relative* conformer
  energies are fairly consistent between implementations. Charges produced by the
  toolkit may generate their own conformers for AM1-BCC; NAGL uses the molecular
  graph. Preserve one charge vector when comparing conformer energies. No universal
  error bound across charge methods is implied.
- `working_mol.to_file(path)` requires an explicit `file_format`; calling it without
  one raises `TypeError: FrozenMolecule.to_file() missing 1 required positional
  argument: 'file_format'`.

## Failure modes and safety

- **Mixed species in one SDF**: the equality check raises
  `ValueError("Multiple chemical species loaded")`. Split the file first.
- **Multi-conformer `Molecule` written to SDF**: `Molecule.to_file(..., file_format="sdf")`
  writes only **one** record (the first conformer) — verified. Write each conformer to
  its own file, as above, or loop manually.
- **Missing NAGL backend**: without `openff-nagl`/`openff-nagl-models` (or with a
  toolkit registry that cannot do ML charges), Sage 2.3.0 charge assignment fails; the
  `NAGLCharges` handler needs the `openff-gnn-am1bcc-1.0.0.pt` model file.
- **Constrained force field by accident**: `openff-2.3.0.offxml` constrains all X-H
  bonds; using it here freezes the very coordinates minimization should relax. Always
  use the `openff_unconstrained-*` file for single points and minimization.
- **Filename safety**: output names come from `molecule.name`, which may be empty or
  contain spaces/slashes. Sanitize it (the notebook falls back to `to_hill_formula()`).
- **Do not "fix" the physics**: leave the force field version, the unconstrained
  choice, and the minimizer defaults (`interchange.minimize(engine="openmm")`) alone
  unless the user asks. Report which force field and charge model produced any numbers
  you hand back.

## Sources

Read at `openforcefield/openff-toolkit` `main` @
`ef074d58a2e7b03625aae7cd3aff87a8ee46c38b` (2026-10-16); the three files below are
byte-identical on `stable` @ `3673715c3223e1d0572a41a346fae113e8aff287` (2026-08-12).

- Notebook: <https://github.com/openforcefield/openff-toolkit/blob/ef074d58a2e7b03625aae7cd3aff87a8ee46c38b/examples/conformer_energies/conformer_energies.ipynb>
- CLI script: <https://github.com/openforcefield/openff-toolkit/blob/ef074d58a2e7b03625aae7cd3aff87a8ee46c38b/examples/conformer_energies/conformer_energies.py>
- Directory: <https://github.com/openforcefield/openff-toolkit/tree/main/examples/conformer_energies>
- `examples/conformer_energies/README.md`, `examples/conformer_energies/.gitignore`,
  `examples/README.md`, `devtools/conda-envs/examples.yaml`, `.github/workflows/examples.yml`
- Force field file `openff-2.3.0.offxml` from `openforcefield/openff-forcefields` @
  `be06586b4b8056dac037278a3b717eef65935bce`; `openff_unconstrained-2.3.0.offxml` ships
  alongside it in the same `openforcefields/offxml/` directory
- Toolkit FAQ (unconstrained force fields, deterministic charges):
  <https://docs.openforcefield.org/projects/toolkit/en/stable/faq.html>
- API surfaces used: `Molecule.from_file` / `Molecule.to_file`,
  `ForceField.create_interchange`, `Interchange.minimize` (`openff-toolkit`,
  `openff-interchange`), `get_openmm_energies` (`openff.interchange.drivers.openmm`)

Upstream licenses: MIT (`openff-toolkit`, `openff-interchange`, `openff-forcefields`).

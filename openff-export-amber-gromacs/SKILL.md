---
name: openff-export-amber-gromacs
description: Export a SMIRNOFF-parametrized OpenFF Topology to AMBER (prmtop/inpcrd/in) and GROMACS (top/gro/mdp) input files via Interchange, then cross-validate single-point energies between OpenMM, sander and gmx. Use when asked to convert an OpenFF/OpenMM system into AMBER or GROMACS files, to run the OpenFF "using SMIRNOFF in Amber or GROMACS" example, or to check that a non-OpenMM export reproduces the OpenMM energy.
---

# Export an OpenFF system to AMBER and GROMACS

## Purpose and scope

Take a chemical system described by an OpenFF `Topology` plus a SMIRNOFF `.offxml`
force field and write engine-native input files for AMBER and GROMACS, without
leaving Python. In scope: small-molecule systems, parametrization, file export, and
single-point energy cross-validation. Out of scope: system preparation (protonation,
solvation, ligand placement), equilibration/production MD protocols, and choice of
force field — those are scientific decisions the user must make.

## Important correction to the directory blurb

The OMSF directory entry and the example's own `README.md` say this conversion goes
"through the ParmEd library". **That is stale.** The current example
(`export_with_interchange.ipynb`) uses `openff-interchange`: it calls
`ForceField.create_interchange()` and `Interchange.to_amber()` / `.to_gromacs()`.
ParmEd is no longer involved. Do not write ParmEd-based instructions for this
workflow. (The example README also misspells its own notebook as
`export_with_interchange.ipynyb`.)

If the user has only an arbitrary OpenMM `System`, this example is not a general
reverse converter. Ask for the original OpenFF topology, coordinates and force field;
do not pretend these can always be reconstructed from a serialized system.

## Required inputs and clarifications

Ask the user before proceeding:

- The **structure** to parametrize, and its provenance. A PDB alone has no bond
  orders or formal charges; OpenFF needs those. Either supply a molecule file that
  carries them (SDF/MOL2) or pass the molecules explicitly as `unique_molecules`.
- The **force field**, by exact filename. Never silently substitute a different
  version — the `.offxml` filename *is* the reproducibility record. The example uses
  `openff-2.3.0.offxml` (Sage line).
- Whether the system is **periodic or vacuum**. This matters: `Interchange.to_gromacs`
  needs a box, and GROMACS does not support a proper vacuum simulation.
- Which **engines** to validate against, and whether AMBER/GROMACS are installed
  (scientific engines are often licence- or size-gated downloads).

## Environment setup and version caveats

Upstream CI installs the examples environment from
`devtools/conda-envs/examples.yaml` in the toolkit repository. A minimal equivalent:

```bash
micromamba create -n openff-export -c conda-forge \
  python=3.12 \
  "openff-toolkit-base=0.19.0" "openff-interchange-base>=0.5" \
  "openff-forcefields>=2025" \
  "openff-amber-ff-ports>=0.0.3" \
  "openff-nagl-models>=2025" \
  "openff-nagl-base>=0.5.3" pytorch-cpu \
  "rdkit=2024" \
  "openmm>=7.6"
micromamba activate openff-export
```

Caveats that bite in practice:

- **Sage 2.3.0 does not use AM1-BCC.** `openff-2.3.0.offxml` declares
  `<NAGLCharges model_file="openff-gnn-am1bcc-1.0.0.pt">`. It therefore needs
  `openff-nagl-models` and a CPU PyTorch install. If those are missing you get a
  charge-assignment failure, not a warning. `openff-2.2.1.offxml` and earlier still
  use `<ToolkitAM1BCC>`.
- The toolkit requires Python >= 3.12 at the `stable` revision read here.
- The notebook's markdown says partial charges are AM1-BCC; that text predates the
  Sage 2.3.0 switch. Trust the `.offxml` file, not the prose.
- AMBER (`sander`, from AmberTools) and GROMACS (`gmx`) are only needed for the
  *validation* step. Export itself is pure Python and in-memory.

## Steps

For the bundled example, download its input into a new working directory (or use
an existing, user-supplied PDB after changing the path in step 1):

```bash
mkdir openff-export-run && cd openff-export-run
curl --fail --location --output 1_cyclohexane_1_ethanol.pdb \
  https://raw.githubusercontent.com/openforcefield/openff-toolkit/3673715c3223e1d0572a41a346fae113e8aff287/examples/using_smirnoff_in_amber_or_gromacs/1_cyclohexane_1_ethanol.pdb
```

Work in this scratch directory. The example ships `1_cyclohexane_1_ethanol.pdb` (one
cyclohexane + one ethanol, 27 atoms, with a `CRYST1` line, so the resulting topology
is periodic — despite the notebook describing it as "in vacuum").

1. Build a `Topology` and load the force field. Pass the molecules as
   `unique_molecules` so bond orders and formal charges are known:

```python
from openff.toolkit import ForceField, Molecule, Topology

ethanol = Molecule.from_smiles("CCO")
cyclohexane = Molecule.from_smiles("C1CCCCC1")

topology = Topology.from_pdb(
    "1_cyclohexane_1_ethanol.pdb",
    unique_molecules=[ethanol, cyclohexane],
)

forcefield = ForceField("openff-2.3.0.offxml")
```

2. Parametrize into an `Interchange`. This is an in-memory conversion; it does not
   shell out to any engine:

```python
interchange = forcefield.create_interchange(topology)
print(repr(interchange))
# Interchange with 7 collections, periodic topology with 27 atoms.
```

3. Write AMBER and GROMACS files:

```python
interchange.to_amber(prefix="system")
interchange.to_gromacs(prefix="system")
```

`to_amber(prefix="system")` writes `system.prmtop`, `system.inpcrd` and
`system_pointenergy.in` (a sander single-point input). `to_gromacs(prefix="system")`
writes `system.top`, `system.gro` and `system_pointenergy.mdp`. Note the
`to_amber` docstring in Interchange 0.5.x wrongly says it produces "foo.top,
foo.gro, and foo_pointenergy.in" — the real names are as above.

Individual writers (`to_prmtop`, `to_inpcrd`, `to_sander_input`, `to_top`, `to_gro`,
`to_gromacs_mdp`) are available if you want to control paths or GROMACS options
(`monolithic`, `hydrogen_mass`, `decimal`).

## Validation and expected outputs

Do not treat "files were written" as success. Compare single-point energies across
engines. The drivers require each engine's executable on `PATH`.

```python
from pprint import pprint
from shutil import which

from openff.interchange.drivers import get_amber_energies, get_gromacs_energies
from openff.interchange.drivers import get_openmm_energies
from openff.interchange.drivers.all import get_summary_data

pprint(get_openmm_energies(interchange).energies)

if which("sander"):
    pprint(get_amber_energies(interchange).energies)
if which("gmx"):
    pprint(get_gromacs_energies(interchange).energies)

# Optional overview: invokes the engines it can find and may need pandas.
# print(get_summary_data(interchange))
```

The upstream notebook stores the following example outputs (kJ/mol). They were not
recomputed for this skill and may reflect an older environment; use them for orientation,
not as golden values for a fresh run:

| Engine  | Bond     | Angle    | Torsion  | vdW      | Electrostatics | RBTorsion |
| ------- | -------- | -------- | -------- | -------- | -------------- | --------- |
| OpenMM  | 0.172642 | 2.738725 | 17.694555 | 11.286877 | -15.070121    | NaN       |
| Amber   | 0.172799 | 2.738846 | 17.694554 | 11.295126 | -15.071186    | NaN       |
| GROMACS | 0.172641 | 2.738762 | 17.694551 | 11.296365 | -15.077801    | 0.0       |

For automation, `EnergyReport.compare(other, tolerances={...})` raises when terms
exceed per-term tolerances (default 1e-3 kJ/mol). The stored example already exceeds
that default for some nonbonded terms. Investigate switching, cutoffs, electrostatics,
and coordinate precision before defining scientifically justified tolerances. Do not
claim strict equivalence or loosen tolerances merely to obtain a passing check.

## Failure modes and safety

- **`SwitchingFunctionNotImplementedWarning`**: the force field requests an 8.0 Å
  switching distance; AMBER has no switching function, so Interchange falls back to a
  hard cut-off and non-bonded interactions differ. This is expected here and is part
  of why the AMBER energies differ slightly. Report it; do not "fix" it by editing the
  force field.
- **`which("sander")` / `which("gmx")` is false**: the driver raises rather than
  skipping. Guard with `which(...)`, or install AmberTools / GROMACS.
- **Vacuum export to GROMACS**: GROMACS needs a box. The related protein example
  assigns a 5 nm pseudo-vacuum box for this reason; a true vacuum system will not
  round-trip.
- **Cut-off / non-bonded mismatch** on `Interchange.combine()` raises
  `CutoffMismatchError` or `SwitchingFunctionMismatchError` — do not force it.
- **Never overwrite the example's input PDB**, and never hand-edit a generated
  `.prmtop`/`.top` and then claim it reproduces OpenMM.
- A matching single-point energy proves the *export layer* is consistent. It says
  nothing about whether the force field, cut-offs, or system preparation are
  appropriate for the scientific question.

## Sources

- Example directory (moving `stable` pointer):
  <https://github.com/openforcefield/openff-toolkit/tree/stable/examples/using_smirnoff_in_amber_or_gromacs>
- Notebook, read at toolkit `stable` = tag `0.19.0`, commit
  `3673715c3223e1d0572a41a346fae113e8aff287` (2026-08-12):
  <https://github.com/openforcefield/openff-toolkit/blob/3673715c3223e1d0572a41a346fae113e8aff287/examples/using_smirnoff_in_amber_or_gromacs/export_with_interchange.ipynb>
- Example README (contains the stale ParmEd wording):
  <https://github.com/openforcefield/openff-toolkit/blob/3673715c3223e1d0572a41a346fae113e8aff287/examples/using_smirnoff_in_amber_or_gromacs/README.md>
- Examples environment file:
  <https://github.com/openforcefield/openff-toolkit/blob/3673715c3223e1d0572a41a346fae113e8aff287/devtools/conda-envs/examples.yaml>
- Interchange export API docs:
  <https://docs.openforcefield.org/projects/interchange/en/stable/using/output.html>
- Interchange source, read at `main` commit
  `389e949317f7d8a35c59699e3d949eb7aa1cffa3`: `openff/interchange/components/interchange.py`
  (`to_amber`, `to_gromacs`, `to_openmm_simulation`) and
  `openff/interchange/drivers/`.
- Force field definitions, read at `openff-forcefields` `main` commit
  `ee065ef1715ba94953f7736e665dd13bceb1e985`:
  `openforcefields/offxml/openff-2.3.0.offxml` (`<NAGLCharges>`),
  `openforcefields/offxml/openff-2.2.1.offxml` (`<ToolkitAM1BCC>`).
- Upstream license: MIT (openff-toolkit, openff-interchange). Attribution only; this
  skill is not a redistribution of upstream files.

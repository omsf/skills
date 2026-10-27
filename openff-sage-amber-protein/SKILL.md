---
name: openff-sage-amber-protein
description: Combine a Sage/SMIRNOFF-parametrized ligand with a protein parametrized by the SMIRNOFF port of Amber ff14SB, using Interchange.combine, then export to OpenMM, AMBER or GROMACS. Use when asked to run the OpenFF "Sage ligand with an Amber protein" / BRD4 benchmark example, to mix an OpenFF small-molecule force field with a protein force field, or to build a protein-ligand complex from separate receptor and ligand files.
---

# Integrating a Sage ligand with an Amber protein

## Purpose and scope

Build a single parametrized system from two independently parametrized pieces: a small
molecule with an OpenFF "Sage" force field, and a protein with the SMIRNOFF port of
Amber ff14SB. The mechanism is `Interchange.combine()`. In scope: parametrizing each
component, combining, and exporting to OpenMM/AMBER/GROMACS. Out of scope: docking,
solvation, equilibration, and free-energy protocols.

## The current example is BRD4, not toluene/T4 lysozyme

The example directory's `README.md` describes
`toluene_in_T4_lysozyme.ipynb` ("Toluene in complex with T4 lysozyme L99A in
TIP3P-FB water"). **That notebook no longer exists** in `stable` or `main`. The only
notebook in the directory is `BRD4_inhibitor_benchmark.ipynb`, which parametrizes
BRD4 inhibitors from the Mobley & Gilson benchmark set. The rendered documentation
also points at `BRD4_inhibitor_benchmark.html`. Follow the notebook, and tell the user
that the README is stale rather than inventing a toluene workflow.

## Required inputs and clarifications

- A **ligand** file that carries bond orders and formal charges (SDF/MOL2) plus a
  matching **receptor** PDB. Both must be in the *same coordinate frame* — this
  workflow combines parameters, it does not place or dock the ligand.
- The **ligand force field** by exact filename (example: `openff-2.3.0.offxml`).
- The **protein force field** by exact filename (example:
  `ff14sb_off_impropers_0.0.4.offxml`).
- Whether the user wants a real solvent box or just a parametrization check. The
  upstream notebook does **not** solvate; it sets a 5 nm "pseudo-vacuum" box purely so
  GROMACS will accept the export.
- Whether this is a benchmark reproduction or a production system. It is a
  parametrization benchmark. Say so explicitly.

## Environment setup and version caveats

```bash
micromamba create -n openff-sage-amber -c conda-forge \
  python=3.12 \
  "openff-toolkit-base=0.19.0" "openff-interchange-base>=0.5" \
  "openff-forcefields>=2025" \
  "openff-amber-ff-ports>=0.0.3" \
  "openff-nagl-models>=2025" \
  "openff-nagl-base>=0.5.3" pytorch-cpu \
  "rdkit=2024" \
  "openmm>=7.6" \
  requests
micromamba activate openff-sage-amber
```

- `ff14sb_off_impropers_0.0.4.offxml` comes from the **`openff-amber-ff-ports`**
  package, not `openff-forcefields`. It has been a dependency of the toolkit since
  0.11.0, so it is usually already present. List what is available with
  `openff.toolkit.typing.engines.smirnoff.forcefield.get_available_force_fields()`.
- Sage 2.3.0 uses `<NAGLCharges>` (a GNN charge model) rather than AM1-BCC, so
  `openff-nagl-models` and CPU PyTorch are required. See the
  `openff-gnn-am1bcc-1.0.0.pt` model. Do not fall back to another charge method.
- The SMIRNOFF ff14SB port is a **proof-of-concept**. Upstream warns it "may not
  precisely match the energetics or forces of the original ff14SB, and it may not be
  appropriate in combination with Sage". It is not production-ready.

## Steps

1. Work in a new scratch directory so input files cannot be overwritten.
   Fetch the benchmark inputs (the URLs below pin the revision read for this skill). The notebook downloads from the MobleyLab
   `benchmarksets` repository:

```python
# https://github.com/MobleyLab/benchmarksets
import requests

repo_url = (
    "https://raw.githubusercontent.com/MobleyLab/benchmarksets/"
    "5233e0a1529b3c3ba9adf8c245de2b658e5fbf40/input_files/"
)
sources = {
    "receptor.pdb": repo_url + "BRD4/pdb/BRD4.pdb",
    "ligand.pdb": repo_url + "BRD4/pdb/ligand-1.pdb",
    "ligand.sdf": repo_url + "BRD4/sdf/ligand-1.sdf",
}
for filename, url in sources.items():
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    with open(filename, "x") as handle:
        handle.write(r.text)
```

   `ligand.pdb` is downloaded but unused by the notebook. `BRD4.pdb` is the apo
   receptor (2035 atoms, no `HETATM`, no `CRYST1`, uses `HIE` for 17 residues);
   `ligand-1.sdf` is the docked pose (26 atoms: 12 C, 9 H, 4 N, 1 S).

2. Parametrize the ligand and the receptor separately, then combine:

```python
from openff.toolkit import ForceField, Molecule, Quantity, Topology

ligand_molecule = Molecule.from_file("ligand.sdf")

sage = ForceField("openff-2.3.0.offxml")
ligand_system = sage.create_interchange(topology=ligand_molecule.to_topology())

receptor_topology = Topology.from_pdb("receptor.pdb")
ff14sb = ForceField("ff14sb_off_impropers_0.0.4.offxml")
receptor_system = ff14sb.create_interchange(topology=receptor_topology)

complex_system = receptor_system.combine(ligand_system)

# Make pseudo-vacuum because GROMACS does not support proper vacuum simulations
complex_system.box = Quantity([5, 5, 5], "nanometer")
```

   `combine()` appends the `other` topology to `self` without reordering atoms, so the
   receptor keeps atom indices `0..2034` and the ligand follows.

3. Export:

```python
complex_system.to_openmm()
complex_system.to_amber("complex")
complex_system.to_gromacs(prefix="complex")
```

   Files: `complex.prmtop`, `complex.inpcrd`, `complex_pointenergy.in`,
   `complex.top`, `complex.gro`, `complex_pointenergy.mdp`.

## Validation and expected outputs

- Confirm the combined atom count is receptor + ligand (2035 + 26 = 2061) and that
  `complex_system.topology.n_atoms` matches before writing anything.
- Confirm the ligand retains its Sage parameters: inspect
  `complex_system["vdW"]` / `complex_system["Electrostatics"]` parameter collections,
  or compare `ligand_system` energies with `get_openmm_energies(ligand_system)`.
- Run the energy drivers on the combined system if the engines are available
  (`openff.interchange.drivers.get_openmm_energies`, `get_amber_energies`,
  `get_gromacs_energies`, or `get_summary_data`). Expect OpenMM, sander and gmx
  single-point energies to be compared term by term in matching units. Switching
  differences and coordinate precision can exceed 1e-3 kJ/mol; investigate and report
  them rather than asserting that all engine settings are identical.
- Two warnings are **expected** and should be reported rather than suppressed:

  - `InterchangeCombinationWarning: Interchange object combination is complex and may
    produce strange results outside of use cases it has been tested in.`
  - `InterchangeCombinationWarning: Found electrostatics 1-4 scaling factors of 5/6
    with slightly different rounding (0.833333 and 0.8333333333).` The value
    `0.8333333333` is used, "which may or may not introduce small errors".

  Also expect `UserWarning: Ambiguous failure while processing constraints.
  Constraining h-bonds as a stopgap.` on export, and
  `SwitchingFunctionNotImplementedWarning` for AMBER.

## Failure modes and safety

- **`CutoffMismatchError` / `SwitchingFunctionMismatchError` / `UnsupportedCombinationError`**
  from `combine()`: the two force fields disagree on non-bonded treatment. Do not
  patch the cut-offs by hand to force the combination through; report the mismatch and
  ask the user how to resolve it.
- **`Topology.from_pdb` on the receptor fails** if atoms are missing or hydrogens are
  not explicit, or if non-standard residue/atom names are used. `BRD4.pdb` is already
  complete. If the user supplies their own receptor, they must prepare it first
  (protonation states, termini, missing loops) — that is their scientific decision.
- **Silently changing force field versions** invalidates the benchmark. Pin the
  filenames and record them in the methods text.
- **A 5 nm box is not solvation.** Simulations in it are pseudo-vacuum and are not
  comparable to solvated free-energy benchmarks. Do not present them as such.
- Ask before any long or GPU-accelerated run. This workflow itself is cheap
  (parametrization only) and safe to run on CPU.

## Sources

- Example directory (moving `stable` pointer):
  <https://github.com/openforcefield/openff-toolkit/tree/stable/examples/using_smirnoff_with_amber_protein_forcefield>
- Notebook, read at toolkit `stable` = tag `0.19.0`, commit
  `3673715c3223e1d0572a41a346fae113e8aff287` (2026-08-12):
  <https://github.com/openforcefield/openff-toolkit/blob/3673715c3223e1d0572a41a346fae113e8aff287/examples/using_smirnoff_with_amber_protein_forcefield/BRD4_inhibitor_benchmark.ipynb>
- Directory README (stale; names a removed notebook):
  <https://github.com/openforcefield/openff-toolkit/blob/3673715c3223e1d0572a41a346fae113e8aff287/examples/using_smirnoff_with_amber_protein_forcefield/README.md>
- Rendered docs page for the same notebook:
  <https://docs.openforcefield.org/en/latest/examples/openforcefield/openff-toolkit/using_smirnoff_with_amber_protein_forcefield/BRD4_inhibitor_benchmark.html>
- `Interchange.combine` implementation, read at Interchange `main` commit
  `389e949317f7d8a35c59699e3d949eb7aa1cffa3`:
  `openff/interchange/components/interchange.py` and
  `openff/interchange/operations/_combine.py`.
- Benchmark inputs: <https://github.com/MobleyLab/benchmarksets/tree/master/input_files/BRD4>
  (repository HEAD `5233e0a1529b3c3ba9adf8c245de2b658e5fbf40` at time of reading; the
  notebook uses a moving `master` branch URL; the snippet above pins this commit).
- Force field ports, read at `openff-amber-ff-ports` `main` commit
  `c2a0c45c92240b6f695059291abaae8df5343913`:
  `openff/amber_ff_ports/offxml/ff14sb_off_impropers_0.0.4.offxml`.
- Upstream licenses: MIT (openff-toolkit, openff-interchange, openff-amber-ff-ports);
  MobleyLab `benchmarksets` data — check the repository for dataset terms before
  redistribution. Attribution only; this skill is not a redistribution of upstream files.

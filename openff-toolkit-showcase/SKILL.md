---
name: openff-toolkit-showcase
description: Prepare, parametrize, simulate and visualize a solvated protein-ligand complex end-to-end with the OpenFF Toolkit, Interchange and OpenMM (the OpenFF Toolkit Showcase). Use when asked to run the toolkit showcase notebook, to go from a prepared receptor PDB plus ligand SDF to a short OpenMM trajectory, or to see the full Topology -> ForceField -> Interchange -> Simulation pipeline in one place.
---

# OpenFF Toolkit Showcase: prepare and run a protein-ligand simulation

## Purpose and scope

The end-to-end tour of the OpenFF stack: assemble a solvated protein-ligand complex,
parametrize it by combining a small-molecule force field (Sage) with a SMIRNOFF port of
Amber ff14SB, convert it to an OpenMM `Simulation` via Interchange, minimize, run a
short trajectory, and visualize it. It is a teaching example, **not** a validated
simulation protocol.

## Required inputs and clarifications

- A **prepared receptor PDB**: hydrogens explicit, termini capped, missing atoms
  rebuilt, crystallographic waters present. The example ships `5tbm_prepared.pdb`
  (1770 atoms, no `CRYST1`) and `PT2385.sdf` (38 atoms, a HIF2α antagonist from the
  Merck `fep-benchmark` set).
- A **ligand SDF** with bond orders, formal charges and coordinates super-imposable on
  the receptor.
- **Force field choice**: the example uses `ForceField("openff-2.3.0.offxml",
  "ff14sb_off_impropers_0.0.4.offxml")`. Order matters — SMIRNOFF applies the *last*
  matching parameter, so the general force field goes first.
- **Simulation intent**: duration, temperature/pressure, and acknowledgement that the
  example has no equilibration. Ask before any run that costs real compute time or
  uses a GPU.

## Environment setup and version caveats

Upstream CI runs this notebook with `devtools/conda-envs/examples.yaml`:

```bash
micromamba create -n openff-showcase -c conda-forge python=3.12 \
  "openff-toolkit-base=0.19.0" "openff-interchange-base>=0.5" "openff-forcefields>=2025" \
  "openff-amber-ff-ports>=0.0.3" "openff-nagl-models>=2025" "rdkit=2024" \
  "openmm>=7.6" mdtraj nglview pdbfixer numpy
micromamba activate openff-showcase
```

Caveats:

- **Not production-ready.** Upstream marks the SMIRNOFF ff14SB port as a
  proof-of-concept: "It may not precisely match the energetics or forces of the
  original ff14SB, and it may not be appropriate in combination with Sage." The
  notebook warns to wait for protein support in the Rosemary force field. Rosemary is
  now available in alpha (see `openff-ptm-prototype`) — switching to it is a
  scientific decision, so ask rather than substituting it silently.
- The notebook's prose still calls Sage "2.0.0" and says charges are AM1-BCC, but the
  code uses `openff-2.3.0.offxml`, which declares `<NAGLCharges>` and therefore needs
  `openff-nagl-models` plus CPU PyTorch. Trust the `.offxml`.
- `sim.mdp` and `minim.mdp` sit in the example directory but are **unused** by the
  current notebook — leftovers from a GROMACS variant. Do not wire them in silently.
- `nglview` is only needed for visualization; the numerics work headless. The notebook's
  `%env OMP_NUM_THREADS=2` is a notebook magic, not Python — set the variable in the
  environment when running as a script.

## Steps

Download the two bundled inputs into a fresh working directory:

```bash
mkdir openff-showcase-run && cd openff-showcase-run
base=https://raw.githubusercontent.com/openforcefield/openff-toolkit/3673715c3223e1d0572a41a346fae113e8aff287/examples/toolkit_showcase
curl --fail --location "$base/5tbm_prepared.pdb" --output 5tbm_prepared.pdb
curl --fail --location "$base/PT2385.sdf" --output PT2385.sdf
export OMP_NUM_THREADS=2
```

1. Load the ligand and the receptor, then solvate the receptor with PDBFixer:

```python
import openmm
import openmm.app
import openmm.unit as openmm_unit
from pdbfixer import PDBFixer
from openff.toolkit import ForceField, Molecule, Topology

receptor_path = "5tbm_prepared.pdb"
ligand_path = "PT2385.sdf"

ligand = Molecule.from_file(ligand_path)

fixer = PDBFixer(receptor_path)
fixer.addSolvent(
    padding=1.0 * openmm_unit.nanometer,
    ionicStrength=0.15 * openmm_unit.molar,
)
with open("receptor_solvated.pdb", "w") as f:
    openmm.app.PDBFile.writeFile(fixer.topology, fixer.positions, f)

top = Topology.from_pdb("receptor_solvated.pdb")
```

   The notebook gets `openmm.app` as a side effect of importing `pdbfixer`; import it
   explicitly as above so the code does not depend on import order.

2. Append the ligand and drop any solvent that clashes with it. This is a compacted
   adaptation of the notebook helper: it refuses to delete protein/cofactors, whereas
   the original can drop any clashing molecule. The 1.5 Å radius is unchanged:

```python
import numpy as np
from openff.units import Quantity

def insert_molecule_and_remove_clashes(
    topology: Topology,
    insert: Molecule,
    radius: Quantity = Quantity(1.5, "angstrom"),
) -> Topology:
    """Append a ligand, removing clashing water or monatomic ions only.

    Parameters
    ----------
    topology : Topology
        Solvated receptor with positions and box vectors.
    insert : Molecule
        Ligand in the receptor coordinate frame, with a conformer.
    radius : Quantity
        Minimum allowed atom distance.

    Returns
    -------
    Topology
        Combined topology with unchanged box vectors.
    """
    water = Molecule.from_smiles("O")
    insert_xyz = insert.conformers[0][:, None, :].m_as("nanometer")
    kept: list[Molecule] = []
    for molecule in topology.molecules:
        xyz = molecule.conformers[0][None, :, :].m_as("nanometer")
        if np.linalg.norm(xyz - insert_xyz, axis=-1).min() > radius.m_as("nanometer"):
            kept.append(molecule)
        else:
            if molecule.n_atoms != 1 and not molecule.is_isomorphic_with(water):
                raise ValueError("Ligand clashes with protein/cofactor; check the pose")
            print(f"Removed {molecule.to_smiles()} molecule")

    kept.append(insert)
    new_topology = Topology.from_molecules(kept)
    new_topology.box_vectors = topology.box_vectors
    return new_topology

top = insert_molecule_and_remove_clashes(top, ligand)
with open("topology.json", "w") as f:
    print(top.to_json(), file=f)
```

   `topology.json` is the durable artifact: molecules, bond orders, formal charges,
   conformers and box vectors, reloadable with `Topology.from_json()`. Reuse it
   instead of re-running solvation.

3. Parametrize and build the simulation:

```python
sage_ff14sb = ForceField("openff-2.3.0.offxml", "ff14sb_off_impropers_0.0.4.offxml")
interchange = sage_ff14sb.create_interchange(top)
```

   Parameterization time depends on the system and environment. Inspect warnings
   before creating the simulation; no timing or convergence is guaranteed.

4. Configure OpenMM, minimize, and run:

```python
integrator = openmm.LangevinIntegrator(
    300 * openmm_unit.kelvin, 1 / openmm_unit.picosecond, 0.002 * openmm_unit.picoseconds
)
simulation = interchange.to_openmm_simulation(
    integrator=integrator, platform=openmm.Platform.getPlatformByName("CPU")
)
stride = 100

simulation.reporters.append(
    openmm.app.DCDReporter(file="trajectory.dcd", reportInterval=stride)
)

simulation.minimizeEnergy()
simulation.context.setVelocitiesToTemperature(300 * openmm_unit.kelvin)
simulation.runForClockTime(1.0 * openmm_unit.minute)
```

   This snippet explicitly selects CPU and a 100-step reporting interval; upstream
   can select a GPU automatically and uses stride 1000 when CUDA is available. Record
   the selected platform and approve GPU use before changing it.
   `to_openmm_simulation` reuses the topology's positions and box vectors. Pass any
   barostat or extra force via `additional_forces=[...]`; adding a force to the `System`
   after context creation requires re-initialization.

5. Visualize (optional, notebook only):

```python
import mdtraj
import nglview

trajectory = mdtraj.load(
    "trajectory.dcd", top=mdtraj.Topology.from_openmm(interchange.to_openmm_topology())
)
view = nglview.show_mdtraj(trajectory.image_molecules())
view.add_representation("line", selection="protein")
```

## Validation and expected outputs

- Print the ligand SMILES before and after assembly; `Molecule.from_file` must not
  silently re-perceive the ligand.
- Check the clash removal: the helper prints one line per removed molecule. It should
  remove only waters/ions, never protein. Any removed protein means the radius or the
  ligand placement is wrong — stop.
- Confirm `top.box_vectors` is not `None` and that
  `interchange.topology.n_atoms == top.n_atoms`.
- Report potential energy (kJ/mol) and maximum per-atom force (kJ/(mol·nm)) before
  and after minimization, as the notebook's `describe_state` helper does.
  Minimization exists to remove clashes and force-field/crystal disagreements; if the
  maximum force does not drop sharply, the system is not ready.
- `trajectory.dcd` should exist and grow; check `len(trajectory)` and that backbone
  RMSD does not explode. A 1-minute trajectory is a smoke test, nothing more.

## Failure modes and safety

- **No equilibration.** The notebook jumps from minimization straight to dynamics. Say
  so out loud; do not present the trajectory as equilibrated or production-ready.
- **An unset `OMP_NUM_THREADS`** lets OpenMM oversubscribe cores on shared machines.
- **PDBFixer adds solvent, not protein hydrogens**, so the input must already have
  them. `Topology.from_pdb` fails on missing atoms or non-standard residue/atom names,
  and does not infer bond orders — it raises on chemistry it cannot identify. For
  modified residues use Pablo instead (see `openff-ptm-prototype`).
- **Do not edit `5tbm_prepared.pdb` or `PT2385.sdf` in place.** Write outputs to a
  scratch directory; the example's `.gitignore` already excludes `receptor.pdb`,
  `receptor_solvated.pdb`, `topology.json`, `*.dcd` and friends.
- Ask the user before extending the run time, adding a barostat/equilibration protocol,
  or using a GPU. The upstream run is short but not CPU-only; OpenMM can
  automatically select an available GPU.

## Sources

- Example directory (moving `stable` pointer):
  <https://github.com/openforcefield/openff-toolkit/tree/stable/examples/toolkit_showcase>
- Notebook, read at toolkit `stable` = tag `0.19.0`, commit `3673715c3223e1d0572a41a346fae113e8aff287`
  (2026-08-12):
  <https://github.com/openforcefield/openff-toolkit/blob/3673715c3223e1d0572a41a346fae113e8aff287/examples/toolkit_showcase/toolkit_showcase.ipynb>
- Also at that commit: `examples/toolkit_showcase/README.md` and
  `devtools/conda-envs/examples.yaml`.
- Upstream data provenance: <https://github.com/MCompChem/fep-benchmark> (`hif2a`
  system; `5tbm_prepared.pdb` and the last ligand of `ligands.sdf` were copied out of
  it, per the notebook's own shell snippet).
- Interchange `to_openmm_simulation`, read at Interchange `main` commit
  `389e949317f7d8a35c59699e3d949eb7aa1cffa3`.
- Upstream licenses: MIT (openff-toolkit, openff-interchange). Attribution only; this
  skill is not a redistribution of upstream files.

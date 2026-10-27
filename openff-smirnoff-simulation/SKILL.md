---
name: openff-smirnoff-simulation
description: Build a SMIRNOFF-parameterized system with the OpenFF Toolkit and run a short OpenMM Langevin dynamics simulation via Interchange, for the openff-toolkit SMIRNOFF_simulation example. Use for requests about parametrizing a PDB/SMILES system with Sage, creating an Interchange, writing OpenMM reporters, or converting an OpenFF system to a simulation engine.
---

# OpenFF Example: SMIRNOFF Simulation

## Purpose and scope

Take a structure (here a small cyclohexane/ethanol box shipped with the toolkit),
assign SMIRNOFF "Sage" parameters, and propagate it with OpenMM's Langevin
integrator through `openff-interchange`. This is a **smoke test / teaching
workflow**, not a production MD protocol: 1000 steps, no barostat, no equilibration
schedule. Scope ends at producing an OpenMM trajectory; see the
`using_smirnoff_in_amber_or_gromacs` example for non-OpenMM engines.

Upstream files: `examples/SMIRNOFF_simulation/run_simulation.ipynb` and
`examples/SMIRNOFF_simulation/README.md`.

## Required inputs and clarifications

1. **System definition**: a PDB path (or your own) plus the concrete `Molecule`
   objects that describe bond orders/charges. PDB files are not a reliable source of
   bond orders, so the notebook supplies molecules from SMILES and passes them as
   `unique_molecules`.
2. **Size**: the example defaults to one ethanol + one cyclohexane (27 atoms). The
   larger 40/60% cyclohexane–ethanol `packmol_boxes` file has 3015 atoms — ask before
   switching, it is orders of magnitude slower.
3. **Force field**: the notebook uses `openff-2.3.0.offxml` (Sage 2.3.0). Do not
   silently change the force field version; different Sage releases have different
   parameters and different default partial-charge models.
4. **Charge model**: Sage 2.3.0 declares `NAGLCharges` with model
   `openff-gnn-am1bcc-1.0.0.pt`; this is a learned charge model, not an on-the-fly
   AM1-BCC calculation. Confirm the exact charge model and any intended charge overrides.
5. **Whether the request is "gas phase"**: the directory summary calls this a gas-phase
   example, but the shipped PDB carries `CRYST1 31.885 31.885 31.885` — a periodic
   cubic box, so `Topology.from_pdb` sets box vectors and OpenMM uses periodic
   nonbonded methods. It is a very dilute box, not a vacuum system.

## Environment setup and version caveats

```bash
mamba create -n openff-toolkit -c conda-forge "openff-toolkit=0.19.0"
mamba activate openff-toolkit
```

- The canonical environment for *running the notebooks* is
  `devtools/conda-envs/examples.yaml` (environment name `openff-toolkit-examples`);
  CI runs `pytest --nbval-lax --dist loadscope --ignore=examples/deprecated examples`.
  It pins `openmm >=7.6`, `openff-forcefields >=2025`, `openff-interchange-base >=0.5`,
  `openff-nagl-base >=0.5.3`, `openff-nagl-models >=2025`.
- `openff-forcefields` supplies the SMIRNOFF/Available-force-field files through the
  `openforcefield.smirnoff_forcefield_directory` entry point. `openff-2.3.0.offxml`
  resolves from that package, not from `openff/toolkit/data/`.
- At this revision the released toolkit is **0.19.0** and `main` is post-0.19.0
  development, so docs on the `stable` docs site may lag `main`.
- RDKit and AmberTools are installed by `openff-toolkit`; `openff-toolkit-base`
  installs only hard dependencies. OpenEye is optional and licensed.
- The notebook kernel metadata references a project-specific env
  (`openff-toolkit-test`); that kernel does not exist for you — select your own.
- The Toolkit provides `Molecule.to_topology()` and `ForceField.create_interchange()`;
  Interchange provides `Interchange.to_openmm*` and OpenMM position conversion.

## Steps

### 1. Load molecules and the box

```python
from openff.toolkit import ForceField, Molecule, Topology
from openff.toolkit.utils import get_data_file_path

# 1 molecule of ethanol and 1 of cyclohexane (27 atoms, 31.885 A cubic box).
pdb_file_path = get_data_file_path("systems/test_systems/1_cyclohexane_1_ethanol.pdb")

# Larger, slower alternative (3015 atoms):
# pdb_file_path = get_data_file_path("systems/packmol_boxes/cyclohexane_ethanol_0.4_0.6.pdb")

ethanol = Molecule.from_smiles("CCO")
cyclohexane = Molecule.from_smiles("C1CCCCC1")
```

Optional OpenEye-only variant: `Molecule.from_file(get_data_file_path("molecules/ethanol.mol2"))`
guarded by `OpenEyeToolkitWrapper.is_available()`.

### 2. Build the Topology and parametrize

```python
topology = Topology.from_pdb(pdb_file_path, unique_molecules=[ethanol, cyclohexane])

forcefield = ForceField("openff-2.3.0.offxml")   # Sage 2.3.0
interchange = forcefield.create_interchange(topology)

assert topology.box_vectors is not None   # CRYST1 was read from the PDB
```

`Topology.from_pdb` warns that molecules adopt the PDB's geometric stereochemistry
even when it conflicts with the stereochemistry in `unique_molecules`; it also
rejects pseudoparticles and requires explicit hydrogens and CONECT records for
non-polymer, non-water molecules.

### 3. Configure and run the simulation

```python
import openmm
import openmm.app
from openmm import unit

time_step = 2 * unit.femtoseconds
temperature = 300 * unit.kelvin
friction = 1 / unit.picosecond
integrator = openmm.LangevinIntegrator(temperature, friction, time_step)

num_steps = 1000
trj_freq = 100
data_freq = 100

simulation = interchange.to_openmm_simulation(
    integrator, platform=openmm.Platform.getPlatformByName("CPU")
)
simulation.context.setVelocitiesToTemperature(temperature)

simulation.reporters.append(openmm.app.PDBReporter("trajectory.pdb", trj_freq))
simulation.reporters.append(
    openmm.app.StateDataReporter(
        "data.csv",
        data_freq,
        step=True,
        potentialEnergy=True,
        temperature=True,
        density=True,
    )
)

simulation.step(num_steps)
```

This snippet explicitly selects CPU rather than letting OpenMM choose a GPU. Record
that platform choice; use another platform only with approval.

### 4. Inspect the outputs

```bash
head data.csv
```

Or in Python (preserve the header; normalize its first field, which starts with `#`):

```python
import pandas as pd

df = pd.read_csv("data.csv")
df = df.rename(columns={df.columns[0]: "Step"})
print(df.columns.tolist())
print(df[df["Step"] % 500 == 0])
```

## Validation and expected outputs

- Files written to the **notebook's working directory**: `trajectory.pdb`
  (a frame every 100 steps) and `data.csv` whose header is exactly
  `#"Step","Potential Energy (kJ/mole)","Temperature (K)","Density (g/mL)"`. The
  stored upstream run reports potential energies of roughly 50–90 kJ/mol.
- `Temperature (K)` is only an instantaneous kinetic estimate for a 27-atom system, so
  it swings widely: the stored upstream run oscillates roughly between 170 K and 360 K
  across frames. Do not treat that spread as a failure; a runaway or NaN temperature
  does indicate bad structures or parameters.
- `Density (g/mL)` is constant here (~0.0067 g/mL for the two-molecule 31.885 Å box) and
  is not a physically meaningful liquid density.
- Round-trip check worth adding: the same `Interchange` can emit an `openmm.System`
  and `openmm.app.Topology`; particle counts in the topology/positions must agree
  (`to_openmm_positions` includes virtual sites by default).

## Failure modes and safety

- **Missing `unique_molecules`**: `Topology.from_pdb` will fail or guess connectivity
  for molecules without CONECT records. PDB-derived bond orders are unreliable by
  design.
- **Stereo mismatch**: geometric stereochemistry from the PDB wins; if that matters for
  your science, re-derive conformers from the `Molecule` objects instead.
- **Charge assignment**: Sage 2.3.0 needs the NAGL runtime, PyTorch, and the named
  model file. Older force fields using `ToolkitAM1BCC` instead need the selected
  semi-empirical backend. Do not substitute charge methods to work around an install error.
- **Interchange not installed / version skew**: `ForceField.create_interchange` imports
  `openff.interchange` lazily and raises if the version does not satisfy the toolkit's
  requirement.
- **Do not "improve" the protocol silently**: 300 K, 2 fs, 1/ps friction, 1000 steps is
  what upstream asserts. Adding a barostat or changing `nonbondedMethod` without
  stating it changes results; extra forces such as a barostat should be supplied via
  `to_openmm_simulation(..., additional_forces=[...])` so the `Context` is built once.
- **No minimization here**: for real systems minimize and equilibrate first. (The
  virtual-sites example is the case where skipping minimization is dangerous.)
- **Resource use**: the 3015-atom box is still small, but a GPU is not needed and not
  requested. Confirm the wall-clock budget for larger systems.

## Sources

Read at commit `ef074d58a2e7b03625aae7cd3aff87a8ee46c38b`
(`openforcefield/openff-toolkit`, `main`, 2026-10-16); the same files are byte-identical
on branch `stable` at `3673715c3223e1d0572a41a346fae113e8aff287` (2026-08-12), so this
skill is valid for both.

- Notebook: `examples/SMIRNOFF_simulation/run_simulation.ipynb`
  <https://github.com/openforcefield/openff-toolkit/blob/ef074d58a2e7b03625aae7cd3aff87a8ee46c38b/examples/SMIRNOFF_simulation/run_simulation.ipynb>
- Directory: <https://github.com/openforcefield/openff-toolkit/tree/main/examples/SMIRNOFF_simulation>
- `examples/SMIRNOFF_simulation/README.md`, `examples/README.md`
- `docs/installation.md`, `docs/users/smirnoff.md`,
  `devtools/conda-envs/examples.yaml`, `.github/workflows/examples.yml`
- Data: `openff/toolkit/data/systems/test_systems/1_cyclohexane_1_ethanol.pdb`,
  `openff/toolkit/data/systems/packmol_boxes/cyclohexane_ethanol_0.4_0.6.pdb`
- API signatures: `ForceField.create_interchange` in
  `openff/toolkit/typing/engines/smirnoff/forcefield.py`,
  `Topology.from_pdb` in `openff/toolkit/topology/topology.py`,
  `Interchange.to_openmm_simulation` in `openff/interchange/components/interchange.py`
  (`openforcefield/openff-interchange` @ `389e949317f7d8a35c59699e3d949eb7aa1cffa3`)
- Force field file: `openff-2.3.0.offxml` from `openforcefield/openff-forcefields`
  @ `ee065ef1715ba94953f7736e665dd13bceb1e985`

Licenses: MIT (`openforcefield/openff-toolkit`, `openff-interchange`,
`openff-forcefields`).

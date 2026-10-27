---
name: openff-ptm-prototype
description: Parametrize and simulate a post-translationally modified protein with OpenFF using the Pablo PDB loader and the Rosemary alpha force field (openff-ptm-prototype). Use when asked to run the PTM prototype notebook, to parameterize a non-canonical residue such as a dye-labelled cysteine, to define a custom ResidueDefinition, or to solvate and simulate a modified protein with graph (NAGL) charges.
---

# Post Translational Modification Prototype

## Purpose and scope

Parameterize a protein containing a **non-canonical residue** (here, a cysteine
labelled with a fluorescein maleimide dye) and run a short simulation. The prototype
combines three things:

1. **Pablo** (`openff-pablo`) — a PDB loader that identifies residues from a
   `ResidueDefinition` library (the CCD by default) instead of guessing.
2. The **Rosemary alpha** force field (`openff_no_water-3.0.0-alpha0.offxml`), which
   covers proteins and small molecules in one SMIRNOFF file.
3. **Graph charges** (NAGL / AshGC), so a large non-canonical residue gets consistent
   partial charges with no library charges.

This is an **explicitly labelled prototype**. Treat it as a demonstration of
infrastructure, not as a validated protocol.

## Required inputs and clarifications

- A **PDB file containing the modified residue**, with all atoms (including
  hydrogens) present and CONECT records for bonds the CCD cannot infer. The repo ships
  `3ip9_dye.pdb` (6184 atoms; 55-atom `DYE` residue at chain A residue 178, 1008
  `HOH`, no `CRYST1` line).
- The **chemistry of the modification**. Pablo needs bond orders and formal charges;
  a PDB does not carry them. Ask the user to describe (or supply) the residue: a
  SMILES, an SDF of the modifying group, or a reaction SMARTS.
- A **force field and water model** by exact filename. The notebook uses
  `ForceField("openff_no_water-3.0.0-alpha0.offxml", "opc3.offxml")`. Do not swap the
  Rosemary alpha for a released Sage force field — Sage has no protein support.
- The intended **simulation**: box shape, padding, salt concentration, ensemble, and
  whether the user accepts an alpha-quality model. Ask before any run beyond a smoke
  test, and before any GPU job.

## Environment setup and version caveats

The repo's `env.yaml` is the authoritative environment. It installs Pablo from GitHub:

```bash
git clone https://github.com/openforcefield/ptm_prototype.git
cd ptm_prototype
git checkout 1db89b53b590dbba5e0940a48b8f9aa0d5788a2e
micromamba env create -n openff_ptm_prototype -c conda-forge -f env.yaml
micromamba install -n openff_ptm_prototype -c conda-forge "packmol>=21" mdtraj
micromamba run -n openff_ptm_prototype jupyter lab ptm_sim.ipynb
```

Caveats, all of which matter:

- **Pablo is pinned to `v0.2.1`** (`pip: git+...openff-pablo.git@v0.2.1`). The README
  claims the repo "has since been updated to install the latest version of Pablo", but
  `env.yaml` at the current `main` commit still pins `v0.2.1`; upstream tag `v0.2.2`
  exists. The APIs used here (`ResidueDefinition.anon_from_sdf`, `.react`,
  `topology_from_pdb`, `STD_CCD_CACHE`) are unchanged between the two, but do not
  "upgrade" silently — record the revision.
- **Rosemary is alpha.** The `openff-forcefields` README lists
  `openff_no_water-3.0.0-alpha0` (2025-10-20) and newer `alpha2b` (2026-09-11). The
  notebook pins `alpha0`. Moving to `alpha2b` changes the model; that is a scientific
  decision, so ask first.
- **`packmol` must be on `PATH`.** `solvate_topology` shells out to it and raises
  `OSError("Packmol not found, cannot run pack_box()")` otherwise. Install the
  `packmol` conda-forge package (>= 21) or get it via AmberTools.
- `openff-nagl-models` plus CPU PyTorch are required: Rosemary alpha uses graph
  charges.
- Pablo's default `STD_CCD_CACHE` does **not** download from the CCD unless
  `STD_CCD_CACHE.auto_download = True`. Canonical amino acids, water and common ions
  are bundled, so this example runs offline. Setting `auto_download = True` reaches
  the wwPDB network — get user consent.
- The notebook imports `solvate_topology` from
  `openff.interchange.components._packmol`, a **private** module. The public home is
  the `openff.packmol` package in newer environments. The reproduction snippet below
  retains the notebook import; migrate only after checking compatible package versions.

## Steps

1. Take a base residue definition and build the modified one by reaction. `react()`
   returns a list of tuples of possible product definitions:

```python
from openff.pablo import STD_CCD_CACHE, ResidueDefinition, topology_from_pdb

cysteine_resdef = STD_CCD_CACHE["CYS"][0]
maleimide_resdef = ResidueDefinition.anon_from_sdf("maleimide.sdf")

# react() takes reactant residue definitions plus reactant/product SMARTS and
# returns every possible outcome as a tuple of product ResidueDefinitions.
dye_resdef = ResidueDefinition.react(
    reactants=[cysteine_resdef, maleimide_resdef],
    reactant_smarts=[
        "[C:10]-[S:1]-[H:2]",
        "[N:3]1-[C:4](=[O:5])-[C:6](-[H:11])=[C:7](-[H:12])-[C:8](=[O:9])-1",
    ],
    product_smarts=[
        "[N:3]1-[C:4](=[O:5])-[C:6](-[H:2])(-[H:11])-[C@:7](-[S:1]-[C:10])(-[H:12])-[C:8](=[O:9])-1",
    ],
)[0][0]
```

   `anon_from_sdf` discards atom names, so the product is matched chemically against
   the PDB rather than by residue/atom name — which is why CONECT records are needed.
   Leaving atoms from the original cysteine are preserved, so the definition describes
   a substructure that can be spliced into a peptide chain.

2. Load the PDB, offering the new definition as a fallback:

```python
topology = topology_from_pdb("3ip9_dye.pdb", additional_definitions=[dye_resdef])
```

3. Solvate:

```python
from openff.interchange.components._packmol import (
    RHOMBIC_DODECAHEDRON, solvate_topology,
)
from openff.units import Quantity

topology = solvate_topology(
    topology,
    nacl_conc=Quantity(0.1, "mol/L"),
    padding=Quantity(1.2, "nm"),
    box_shape=RHOMBIC_DODECAHEDRON,
)
```

   The box is built at a target density with PACKMOL's tolerance, so it is **not**
   equilibrated. `padding` is required here because `3ip9_dye.pdb` has no box vectors.

4. Parametrize:

```python
from openff.toolkit import ForceField

rosemary = ForceField("openff_no_water-3.0.0-alpha0.offxml", "opc3.offxml")
interchange = rosemary.create_interchange(topology)
```

   Expect a `Proton(s) added/removed` warning from NAGL on large molecules. Upstream
   documents it as a harmless, cosmetic bug
   (<https://github.com/openforcefield/openff-nagl/issues/198>).

5. Simulate in OpenMM, serializing the `System` for the record:

```python
import openmm
import openmm.app

simulation = interchange.to_openmm_simulation(
    integrator=openmm.LangevinMiddleIntegrator(
        300 * openmm.unit.kelvin, 1 / openmm.unit.picosecond, 2 * openmm.unit.femtosecond
    ),
    additional_forces=[
        openmm.MonteCarloBarostat(1 * openmm.unit.bar, 300 * openmm.unit.kelvin, 25)
    ],
    platform=openmm.Platform.getPlatformByName("CPU"),
)
simulation.reporters.append(openmm.app.DCDReporter("trajectory.dcd", 1000))

with open("system.xml", "w") as f:
    f.write(openmm.XmlSerializer.serialize(simulation.system))

simulation.context.computeVirtualSites()
simulation.minimizeEnergy()
simulation.context.setVelocitiesToTemperature(simulation.integrator.getTemperature())
simulation.runForClockTime(1.0 * openmm.unit.minute)
```

   CPU is selected explicitly to avoid unapproved GPU use; record this platform.
   `import openmm.app` is needed for `DCDReporter`; the notebook gets it indirectly.
   `computeVirtualSites()` must run before minimization when the force field defines
   virtual sites.

## Validation and expected outputs

- **Residue integration**: visualize the topology and confirm the `DYE` residue is
  bonded into the peptide chain, not floating as a separate ligand, and appears once.
- **Atom bookkeeping**: record `topology.n_atoms` before and after `solvate_topology`;
  the increase must equal waters + ions. `interchange.topology.n_atoms` must match.
- **Charges**: the system must be neutral (or have the intended net charge), and the
  modified-residue atoms in the connected protein must have assigned charges.
  `DYE` is not a separate small-molecule component.
- **Box sanity**: the rhombic dodecahedral box must separate the solute from its
  periodic images by at least twice the non-bonded cut-off.
- **Minimization**: maximum force must drop substantially; a PACKMOL box routinely
  starts with large forces.
- Record `system.xml`, the exact `.offxml` filenames, the Pablo revision, the OpenMM
  platform, and the run length alongside the trajectory.

## Failure modes and safety

- **`topology_from_pdb` cannot match the residue**: usually missing CONECT records,
  missing hydrogens, or a mismatched reaction SMARTS. Read the `PabloError` — it names
  the residue. Do not paper over it by loosening the definition's validation.
- **`KeyError` from `STD_CCD_CACHE[...]`**: the residue is not bundled and
  `auto_download` is falsy. Either enable download (with consent) or supply your own
  definition.
- **`OSError` / `PACKMOLRuntimeError`** from solvation: `packmol` binary missing, or
  packing failed to converge. Try a smaller `tolerance` or a larger box; do not
  shrink the padding below the non-bonded cut-off.
- **The box needs NPT equilibration** before any production run. The notebook says so
  explicitly; repeat it to the user.
- **Do not present Rosemary-alpha results as converged or predictive.** The force
  field is alpha and "hasn't yet been updated to preserve performance on small
  molecules".
- Ask before: enabling CCD downloads, running beyond a ~1 minute smoke test, using a
  GPU, or sharing a private protein structure with any hosted service.
- Never edit `3ip9_dye.pdb` or `maleimide.sdf` in place; generate definitions and
  trajectories in a scratch directory.

## Sources

- Repository (moving `main` pointer): <https://github.com/openforcefield/ptm_prototype>
- Notebook: <https://github.com/openforcefield/ptm_prototype/blob/1db89b53b590dbba5e0940a48b8f9aa0d5788a2e/ptm_sim.ipynb>
- Notebook, environment and README, read at `main` commit
  `1db89b53b590dbba5e0940a48b8f9aa0d5788a2e` (2026-05-07): `ptm_sim.ipynb`,
  `env.yaml` (pins `openff-pablo@v0.2.1`), `README.md`, `maleimide.sdf`,
  `generate_rxn_image.py`, `3ip9_dye.pdb`.
- Pablo: <https://github.com/openforcefield/openff-pablo> — `v0.2.1` commit
  `2c54efd941754c55979c56db7bffb06e804415a0`; `v0.2.2` commit
  `fb364e9a9c0bdbde1262d1e822fe13d823033f40`. APIs verified in
  `openff/pablo/residue.py` (`anon_from_sdf`, `react`, `visualize`) and
  `openff/pablo/_pdb.py` (`topology_from_pdb`). Docs:
  <https://openff-pablo.readthedocs.io/en/stable/>
- `openff.packmol` — `solvate_topology` signature and `packmol` lookup, read at `main`
  commit `3c4049d742a8aa31ae15e8c7cdcd7fb2df6809ed`: `openff/packmol/_packmol.py`.
- Force fields, read at `openff-forcefields` `main` commit
  `ee065ef1715ba94953f7736e665dd13bceb1e985`: `openff_no_water-3.0.0-alpha0.offxml`
  (Rosemary, 2025-10-20), `openff_no_water-3.0.0-alpha2b.offxml` (Rosemary,
  2026-09-11), `opc3.offxml`.
- Interchange `to_openmm_simulation` / `additional_forces`, read at `main` commit
  `389e949317f7d8a35c59699e3d949eb7aa1cffa3`.
- Upstream licenses: MIT (ptm_prototype, openff-pablo, openff-packmol,
  openff-forcefields, openff-interchange). Attribution only; this skill is not a
  redistribution of upstream files.

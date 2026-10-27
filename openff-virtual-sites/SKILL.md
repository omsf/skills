---
name: openff-virtual-sites
description: Specify SMIRNOFF virtual sites (off-atom charges) with the OpenFF Toolkit, apply them through Interchange to an OpenMM system, and validate a SMIRNOFF TIP5P implementation against OpenMM's own TIP5P. Use for requests about virtual sites, lone-pair or bond-charge particles, 4/5/6-point water models, TIP5P, `VirtualSites` OFFXML sections, or the openff-toolkit vsite_showcase example.
---

# Virtual Sites Showcase (SMIRNOFF off-atom charges + TIP5P validation)

## Purpose and scope

Two things: (1) write a SMIRNOFF `VirtualSites` section, append it to Sage, apply it to
a ligand with divalent sulfur and aryl Cl/F substituents, and simulate; (2) numerically
compare a SMIRNOFF implementation of TIP5P water with OpenMM's native TIP5P in a small
water box, asserting that geometry and energy agree. Upstream file:
`examples/virtual_sites/vsite_showcase.ipynb`.

Virtual sites are treated as *features of the potential energy function*, not of the
chemical system. The toolkit distinguishes a virtual **site** (a parameter) from a
virtual **particle** (a massless particle placed by that parameter), and only one site
may place multiple particles.

## Required inputs and clarifications

1. **Chemistry and target motif**: which functional group needs off-atom charge (the
   example uses `[*:2]-[#16X2:1]-[*:3]` for a divalent, two-coordinate sulfur — the
   molecule's thioether — plus `BondCharge` sites on aryl `Cl` and `F`). Get the
   intended SMIRKS and geometry.
2. **Geometry parameter source**: the notebook's numeric choices
   (`distance`, `outOfPlaneAngle`, `charge_increment*`, `sigma`, `epsilon`) come from
   prior force-field development work, not from Sage. Never invent or "tidy" these —
   ask the user for the values and their provenance.
3. **Simulation intent**: a 10 000-step demonstration run is cheap on CPU; anything
   larger (water boxes, solvated ligands) should be confirmed before running.
4. **Whether the user wants both parts or only one.** Part 2 (TIP5P comparison) is
   self-contained and CPU-only.

## Environment setup and version caveats

```bash
mamba create -n openff-toolkit -c conda-forge "openff-toolkit=0.19.0"
mamba activate openff-toolkit
```

- Canonical notebook environment: `devtools/conda-envs/examples.yaml`
  (`openff-toolkit-examples`), which pins `openmm >=7.6`,
  `openff-forcefields >=2025`, `openff-interchange-base >=0.5`, `openff-nagl-*`, plus
  `nglview` and `mdtraj` used only for the optional trajectory visualization.
- The notebook's `interchange.collections["VirtualSites"]` requires the Interchange
  virtual-site support in `openff/interchange/smirnoff` (present in
  `openff-interchange-base`); the SMIRNOFF `VirtualSites` handler is not a plugin.
- `tip5p.offxml` is **not** in `openff/toolkit/data`; it ships with the
  `openff-forcefields` package, resolved through the
  `openforcefield.smirnoff_forcefield_directory` entry point. `openmm.app.ForceField("tip5p.xml")`
  resolves against OpenMM's own `openmm/app/data` directory.
- Toolkit unreleased `main` vs released `0.19.0`: this notebook is on `main` and is
  identical on `stable`, but external API docs may lag.
- Only the `parents` exclusion policy is supported by the toolkit/Interchange; any
  other policy in a `VirtualSites` section raises `UnsupportedExportError` on export.
- `sigma="0.1*angstrom"` etc. in the notebook's custom section are tutorial-scale
  placeholders, deliberately different from the 10 Å used by TIP5P. Do not copy them
  into real work without understanding them.

## Steps

### 1. Define the virtual sites and append to Sage

```python
vsite_offxml = """<?xml version="1.0" encoding="utf-8"?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <VirtualSites version="0.3">
        <VirtualSite
            type="DivalentLonePair"
            name="EP"
            smirks="[*:2]-[#16X2:1]-[*:3]"
            distance="0.70 * angstrom"
            charge_increment1="0.1205*elementary_charge"
            charge_increment2="0.0*elementary_charge"
            charge_increment3="0.1205*elementary_charge"
            sigma="0.1*angstrom"
            epsilon="0.0*kilocalories_per_mole"
            outOfPlaneAngle="54.71384225*degree"
            match="all_permutations" >
        </VirtualSite>
        <VirtualSite
            type="BondCharge"
            name="EP"
            smirks="[*:2][Cl:1]"
            distance="0.4*angstrom"
            charge_increment1="0.2*elementary_charge"
            charge_increment2="0.0*elementary_charge"
            sigma="0.1*angstrom"
            epsilon="0.05*kilocalories_per_mole"
            match="all_permutations" >
        </VirtualSite>
        <VirtualSite
            type="BondCharge"
            name="EP"
            smirks="[*:2][F:1]"
            distance="0.4*angstrom"
            charge_increment1="0.2*elementary_charge"
            charge_increment2="0.0*elementary_charge"
            sigma="0.1*angstrom"
            epsilon="0.05*kilocalories_per_mole"
            match="all_permutations" >
        </VirtualSite>
    </VirtualSites>
</SMIRNOFF>
"""

from openff.toolkit import ForceField

force_field = ForceField("openff-2.3.0.offxml", vsite_offxml)
```

Index `:1` is always the **parent** atom. Parameters with different `name` values can
place several sites on one parent; with the same `name`, the last match wins. Use
`match="all_permutations"` for symmetric matches that must produce multiple particles
(e.g. two lone pairs) and `match="once"` when ordering is ambiguous or a single
particle is wanted.

### 2. Parameterize and count particles

```python
from openff.toolkit import Molecule

molecule = Molecule.from_smiles("c1cc(Cl)ccc1C(=O)CS[C]1=CO[C](F)(F)CC1")
molecule.generate_conformers(n_conformers=1)

interchange = force_field.create_interchange(topology=molecule.to_topology())

assert "VirtualSites" in interchange.collections.keys()
n_virtual_sites = len(interchange.collections["VirtualSites"].key_map)
print(f"There are {n_virtual_sites} virtual particles in this topology.")
```

Sage 2.3.0 assigns NAGL charges using `openff-gnn-am1bcc-1.0.0.pt` before applying
the virtual-site charge redistribution. Ensure its runtime and model are installed.

### 3. Prepare a simulation safely

`minimizeEnergy()` and `computeVirtualSites()` must happen **before**
`setVelocitiesToTemperature()`. Virtual sites started at `[0, 0, 0]` produce huge
initial forces, and high initial virtual-site velocities crash OpenMM (see OpenMM issue
3736 referenced in the notebook).

```python
import openmm
import openmm.app
import openmm.unit
from openff.interchange import Interchange


def prepare_simulation(interchange: Interchange) -> openmm.app.Simulation:
    integrator = openmm.LangevinIntegrator(
        300 * openmm.unit.kelvin,
        1 / openmm.unit.picosecond,
        2 * openmm.unit.femtoseconds,
    )
    simulation = interchange.to_openmm_simulation(
        integrator=integrator, platform=openmm.Platform.getPlatformByName("CPU")
    )

    simulation.minimizeEnergy()
    simulation.context.computeVirtualSites()
    simulation.context.setVelocitiesToTemperature(300 * openmm.unit.kelvin)

    simulation.reporters.append(openmm.app.PDBReporter("trajectory.pdb", 100))
    simulation.reporters.append(
        openmm.app.StateDataReporter(
            "data.csv", 100, step=True, potentialEnergy=True, temperature=True, density=True
        )
    )
    return simulation


simulation = prepare_simulation(interchange)
simulation.step(10000)
```

The snippet selects CPU explicitly; record this platform choice rather than letting
OpenMM silently select a GPU. Afterwards `trajectory.pdb` and `data.csv` are in the working directory. The notebook
visualizes with `mdtraj` + `nglview`; expect spurious bonds drawn between atoms and
virtual sites (virtual sites have no element/type in the PDB).

### 4. TIP5P: SMIRNOFF vs OpenMM (validation)

```python
import openmm.app
from openff.toolkit import ForceField

openff_force_field = ForceField("tip5p.offxml")            # openff-forcefields
openmm_force_field = openmm.app.ForceField("tip5p.xml")    # openmm/app/data
```

Build a lattice of `Molecule.from_mapped_smiles("[O:1]([H:2])[H:3]")` waters (ideal
geometry: 0.9572 Å, 104.52°) as the notebook's `build_water_lattice` does, then compare
`evaluate_openff(...)` against `evaluate_openmm(...)` at `num_duplicates=(2, 2, 2)` and
`spacing=(3.0, 3.0, 3.0)` (Å).

The two code paths differ by design and must be compared in matching particle order:

- OpenMM's TIP5P uses `openmm.OutOfPlaneSite`; the SMIRNOFF/OpenFF path uses
  `openmm.LocalCoordinatesSite` with SMIRNOFF-specified geometry, so small geometric
  differences are expected.
- **Ordering**: OpenMM interleaves virtual particles at the end of *each molecule*;
  OpenFF/Interchange appends all of them to the end of the `System`/`Topology`. The
  notebook's `reorder()` helper remaps indices; `Interchange.to_openmm_topology(collate=True)`
  is the supported alternative to hand-rolled reindexing.
- `to_openmm_positions(interchange)` includes virtual sites by default
  (`include_virtual_sites=True`).

## Validation and expected outputs

Part 1: `"VirtualSites" in interchange.collections` is true, the printed virtual-particle
count matches the number of matched SMIRKS parents (sulfur lone pair + Cl + F), and the
10 000-step run completes with a finite potential energy in `data.csv`.

Part 2: the notebook asserts, per molecule averaged over 8 waters, a coordinate norm
difference < `1e-6 angstrom` and an energy difference
< `1e-4 * openmm.unit.kilojoule_per_mole`, then prints both values. If the assertions
fail, check particle ordering first (the usual cause), then whether both force fields
implement the same TIP5P variant. Exact numerical reproduction across OpenMM versions
is not guaranteed.

## Failure modes and safety

- **Crash on step 0**: virtual sites without minimization/`computeVirtualSites()` —
  see step 3.
- **`UnsupportedExportError`** for exclusion policies other than `parents`.
- **Wrong particle ordering** silently corrupts position/energy comparisons; verify
  counts and ordering explicitly rather than trusting `zip`.
- **Charges**: `charge_increment1`, `charge_increment2`, etc. are increments on the
  corresponding tagged orientation atoms (`:1`, `:2`, ...). The virtual particle gets
  the negative sum of those increments, conserving total charge. They are not
  virtual-particle indices. Check the complete redistributed charge sum.
- **Symmetry double-counting**: `match="all_permutations"` on a symmetric motif can
  place coincident particles (e.g. a 4-point water written as `DivalentLonePair`); use
  `match="once"` or distinct `name`s as appropriate.
- **`Molecule.visualize()` and `Topology.visualize()`** need `nglview` (and conformers
  for molecules); non-interactive/headless runs should skip visualization.
- **Resource guard**: the example runs are short; do not scale the water box or
  ligand simulation into a production run without explicit approval, and never
  reinterpret the tutorial's placeholder `sigma`/`epsilon` values as fitted parameters.

## Sources

Read at commit `ef074d58a2e7b03625aae7cd3aff87a8ee46c38b`
(`openforcefield/openff-toolkit`, `main`, 2026-10-16); the example directory is
byte-identical on branch `stable` at `3673715c3223e1d0572a41a346fae113e8aff287`.

- Notebook: `examples/virtual_sites/vsite_showcase.ipynb`
  <https://github.com/openforcefield/openff-toolkit/blob/ef074d58a2e7b03625aae7cd3aff87a8ee46c38b/examples/virtual_sites/vsite_showcase.ipynb>
- Directory: <https://github.com/openforcefield/openff-toolkit/tree/stable/examples/virtual_sites>
- Docs: `docs/users/virtualsites.md` (handling of `match`, `name`, exclusion policy,
  particle/site distinction, atom ordering)
- SMIRNOFF specification, `VirtualSites` section:
  <https://openforcefield.github.io/standards/standards/smirnoff/#virtualsites-virtual-sites-for-off-atom-charges>
- Force field data: `tip5p.offxml` from `openforcefield/openff-forcefields`
  @ `ee065ef1715ba94953f7736e665dd13bceb1e985`; `tip5p.xml` from
  `openmm/openmm` `wrappers/python/openmm/app/data/tip5p.xml`
- API: `Interchange.to_openmm_simulation`, `to_openmm_topology(collate=...)` in
  `openff/interchange/components/interchange.py`,
  `to_openmm_positions` in `openff/interchange/interop/openmm/_positions.py`,
  `_check_virtual_site_exclusion_policy` in `openff/interchange/interop/common.py`
  (`openforcefield/openff-interchange` @ `389e949317f7d8a35c59699e3d949eb7aa1cffa3`)
- OpenMM virtual-site crash context: <https://github.com/openmm/openmm/issues/3736>

Licenses: MIT (OpenFF Toolkit, Interchange, openff-forcefields); OpenMM is MIT-licensed
scientific software.

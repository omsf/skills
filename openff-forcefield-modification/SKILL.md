---
name: openff-forcefield-modification
description: Edit SMIRNOFF parameters after loading a ForceField with the OpenFF Toolkit — change bond/angle values, add a new SMIRKS-specific parameter, tag cosmetic attributes, and measure the energy effect via OpenMM minimization — for the openff-toolkit forcefield_modification example. Use for force field optimization experiments, adding override parameters, or writing a modified .offxml.
---

# OpenFF Example: Forcefield Modification

## Purpose and scope

Load a Sage force field, label the parameters it assigns to a ligand, then mutate
those parameters in the `ForceField` object and observe the effect on the potential
energy and minimized geometry. Covers: editing an existing parameter, adding a new
more-specific SMIRKS parameter that overrides a generic one, and attaching cosmetic
attributes that survive a write to disk.

This is an **exploration / force-field-fitting** workflow, not a validated force
field: anything produced here is a modified copy of Sage and must be labelled as such.
Scope ends at "modified `.offxml` + before/after energies"; it does not re-fit
parameters against QM data. Upstream files:
`examples/forcefield_modification/forcefield_modification.ipynb` and `README.md`.

## Required inputs and clarifications

1. **Which molecule**: the notebook uses a SMILES-derived ligand
   `CC(C)(C)c1c(O)c(O)c2c(c1O)[C@H]1OCCC[C@H]1[C@H](c1cc(O)c(O)c(F)c1)N2` with one
   generated conformer. Ask whether to substitute the user's molecule and whether
   conformers already exist.
2. **Which parameters and why**: never guess a SMIRKS. Derive the applied parameter
   from `label_molecules()` output (see step 2) and confirm the target with the user
   before changing a scientific value.
3. **Force field**: `openff-2.3.0.offxml` (Sage 2.3.0, constrained). Changing the
   version changes every parameter ID and many SMIRKS.
4. **Do not silently upgrade scientific settings.** The notebook's deliberate
   "nonsense" values (a 10 A C-F bond, a 50 degree H-X-H angle, 179 degree ring angles)
   exist to make the effect visible. Say so explicitly if you carry any into a real
   workflow. Also note the upstream helper calls `molecule.visualize(backend="nglview")`,
   which needs a live Jupyter/nglview session — use the headless variant below.

## Environment setup and version caveats

```bash
mamba create -n openff-toolkit -c conda-forge "openff-toolkit=0.19.0"
mamba activate openff-toolkit
```

- Notebook environment: `devtools/conda-envs/examples.yaml` (adds `nglview`, `qcportal`,
  `qcelemental`, `openff-nagl-models`, `gromacs`, `pdbfixer`).
- `openff.units` provides the `unit` registry and `openff.units.openmm` the
  `to_openmm` / `from_openmm` position converters; do not hand `openmm.unit`
  quantities to the toolkit.
- **Unit trap (verified):** `AngleHandler.AngleType.k` is *declared* as
  `kilocalorie / mole / degree**2`, but Sage's `.offxml` files and the notebook write it
  as `kilocalorie / mole / radian**2`. Pint converts silently, so a number is only
  meaningful with its unit. `134.5019777341 kcal/mol/rad**2` becomes
  `0.0409716454 kcal/mol/deg**2`; the generic Sage `a28` oxygen angle
  (`[*:1]-[#8:2]-[*:3]`) is `239.1104681642 kcal/mol/rad**2`. When copying a number from
  an `.offxml`, copy its unit string too.
- Sage 2.3.0 charges come from `<NAGLCharges>` (needs `openff-nagl` +
  `openff-nagl-models`), so `ForceField.get_partial_charges` is the slow step and is
  worth caching.
- Verified against `openff-toolkit` 0.19.0, `openff-interchange` 0.5.5, `openff-units`
  0.4.0, OpenMM 8.6.1, `openff-2.3.0.offxml` sha256
  `7a0d7195a4b717e29fe14aa2cf1bc3203e87ce4aef9cc2ef10fb5145b7edcf63`.

## Steps

### 1. Build the ligand and cache charges

```python
from openff.units import unit
from openff.units.openmm import from_openmm, to_openmm
from openmm import LangevinIntegrator, Platform, app
from openmm import unit as openmm_unit

from openff.toolkit import ForceField, Molecule

ligand_smiles = "CC(C)(C)c1c(O)c(O)c2c(c1O)[C@H]1OCCC[C@H]1[C@H](c1cc(O)c(O)c(F)c1)N2"
ligand = Molecule.from_smiles(ligand_smiles)
ligand.generate_conformers(n_conformers=1)

force_field = ForceField("openff-2.3.0.offxml")
# Charges are expensive and independent of the parameters we are about to edit
ligand.partial_charges = force_field.get_partial_charges(ligand)
```

### 2. Headless minimize-and-report helper

`charge_from_molecules` is essential for this fixed-charge comparison: it stops every re-parameterization from
recomputing charges.

```python
def minimize_and_report(molecule, forcefield, label):
    topology = molecule.to_topology()
    system = forcefield.create_openmm_system(
        topology, charge_from_molecules=[molecule]
    )
    integrator = LangevinIntegrator(
        300 * openmm_unit.kelvin, 1 / openmm_unit.picosecond,
        0.002 * openmm_unit.picoseconds,
    )
    simulation = app.Simulation(
        topology.to_openmm(), system, integrator,
        platform=Platform.getPlatformByName("CPU"),
    )
    simulation.context.setPositions(to_openmm(molecule.conformers[0]))

    initial = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    simulation.minimizeEnergy()
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    print(f"[{label}] initial {initial.format('%0.1F')}; "
          f"minimized {state.getPotentialEnergy().format('%0.1F')}")
    return from_openmm(state.getPositions(asNumpy=True))

baseline_coordinates = minimize_and_report(ligand, force_field, "baseline")
```

For the notebook's picture instead of the print, add `from copy import deepcopy`, set the
minimized coordinates on a `deepcopy` of the molecule, and call
`vis_mol.visualize(backend="nglview")`.

### 3. Find the parameter you actually want

```python
labels = force_field.label_molecules(ligand.to_topology())

fluorine = next(atom for atom in ligand.atoms if atom.symbol == "F")
carbon = next(iter(fluorine.bonded_atoms))
fc_indices = (fluorine.molecule_atom_index, carbon.molecule_atom_index)
print(labels[0]["Bonds"][fc_indices])   # BondType with .smirks, .length, .k

for atoms, parameter in labels[0]["Angles"].items():   # every H-X-H angle: one param
    symbols = [ligand.atoms[i].symbol for i in atoms]
    if symbols[0] == "H" and symbols[2] == "H":
        print(atoms, parameter.smirks, parameter.id)
```

### 4. Edit an existing parameter

Editing the object returned by `label_molecules` does nothing — edit the copy owned
by the `ForceField`:

```python
fluorine_bond = force_field["Bonds"].parameters["[#6:1]-[#9:2]"]
fluorine_bond.length = 10 * unit.angstrom          # deliberately absurd
minimize_and_report(ligand, force_field, "C-F length edit")

hxh_angle = force_field["Angles"].parameters["[#1:1]-[#6X4:2]-[#1:3]"]
hxh_angle.angle = 50 * unit.degree                 # deliberately absurd
minimize_and_report(ligand, force_field, "C-F and H-X-H edits")
```

`force_field["Angles"].parameters` is a `ParameterList` keyed by SMIRKS; a missing key
raises `KeyError`. Precedence rule: **the last matching parameter in the section wins**,
so generic SMIRKS sit near the top of a section and overrides near the bottom.

### 5. Add a new, more specific parameter

The generic `[*:1]-[#8:2]-[*:3]` matches both ring and hydroxyl oxygens, so appending a
pyran-ring-specific parameter (`add_parameter` appends to the **end** of the section by
default, giving it highest precedence) overrides only the ring oxygens:

```python
import openff.toolkit.typing.engines.smirnoff.parameters as offtk_parameters

# Rebuild from a clean force field so earlier edits do not leak in
force_field = ForceField("openff-2.3.0.offxml")
cyclic_nitrogen_smirks = "[*:1]~[#7X3$(*~[#6X3,#6X2,#7X2+0]):2]~[*:3]"
force_field["Angles"].parameters[cyclic_nitrogen_smirks].angle = 179 * unit.degree

cyclic_oxygen_smirks = "C1C[C:1]-[O:2]-[C:3]C1"
angles_handler = force_field["Angles"]
angles_handler.add_parameter(
    parameter=offtk_parameters.AngleHandler.AngleType(
        smirks=cyclic_oxygen_smirks,
        angle=179 * unit.degree,
        # Sage writes k in radian**-2; copy the unit from the .offxml you are mirroring
        k=134.5019777341 * unit.kilocalorie / (unit.mole * unit.radians**2),
    )
)
```

`AngleType` requires both `angle` and `k`; omitting either raises `SMIRNOFFSpecError`.
Reusing an existing SMIRKS raises `DuplicateParameterError` unless you pass
`allow_duplicate_smirks=True`.

### 6. Cosmetic attributes and writing to disk

```python
angles_handler.parameters[cyclic_nitrogen_smirks].add_cosmetic_attribute(
    attr_name="note", attr_value="funky pseudo-linear cyclic nitrogen"
)
angles_handler.parameters[cyclic_oxygen_smirks].add_cosmetic_attribute(
    attr_name="note", attr_value="funky pseudo-linear cyclic oxygen"
)

force_field.to_file("modified.offxml")
# force_field.to_file("modified.offxml", discard_cosmetic_attributes=True)  # drops note=
```

```bash
grep -o 'note="[^"]*"' modified.offxml
```

**Reload trap (verified):** reading that file back needs
`allow_cosmetic_attributes=True`, otherwise you get
`SMIRNOFFSpecError: Unexpected kwarg (note: ...) passed to ... AngleType constructor`.

```python
reloaded = ForceField("modified.offxml", allow_cosmetic_attributes=True)
print(reloaded["Angles"].parameters[cyclic_oxygen_smirks].k)
print(reloaded["Angles"].parameters[cyclic_oxygen_smirks]._note)
```

## Validation and expected outputs

- Helper output is one `initial`/`minimized` pair per edit, in kJ/mol. Reproduced on
  the notebook ligand (Sage 2.3.0, OpenMM 8.6.1, one generated conformer):
  baseline `initial 361.6 kJ/mol -> minimized -55.2 kJ/mol`;
  after `length = 10 A` on the C-F bond `initial ~111385 -> minimized ~-51.6`;
  after `H-X-H = 50 deg` `initial ~112425 -> minimized ~83.8`.
  These values depend on the generated conformer and minimizer. A parameter can
  match many terms; neither small nor large energy shifts alone establish correctness.
  The *initial* energy is deterministic for a given force field and geometry; the
  minimized value is not (OpenMM minimization stops at a 10 kJ/mol/nm force tolerance),
  so compare minimized numbers only to a few tenths of a kJ/mol.
- `force_field.to_file("modified.offxml")` must contain both `note="..."` attributes;
  with `discard_cosmetic_attributes=True` it must contain zero.
- Reload check: the new `AngleType` round-trips with `k = 134.5019777341
  kilocalorie / mole / radian ** 2`, and `_cosmetic_attribs == ['note']`.
- Sanity check that the override took effect: `label_molecules` on the modified force
  field should now map the ring C-O-C angle to your new SMIRKS and the hydroxyl
  oxygens to the generic one.

## Failure modes and safety

- **Editing the labeled parameter**: `labels[0]["Angles"][idx].angle = ...` has no effect
  on the `ForceField`; always go through `force_field[<tag>].parameters[smirks]`.
- **Generic SMIRKS hitting more atoms than intended**: the notebook's first attempt at
  "ring heteroatoms" also linearizes hydroxyls, because `[*:1]-[#8:2]-[*:3]` matches any
  X-O-X. Print the matched atom indices before committing a change.
- **Stale `System`**: a parameter edit does not invalidate an already-built OpenMM
  `System`; the helper rebuilds it each call, so do the same if you cache one.
- **Mixed charge models**: `charge_from_molecules` reuses whatever charges the molecule
  carries. If they came from a different charge model than the edited force field, you
  silently mix models — recompute after touching `<NAGLCharges>`/`<ToolkitAM1BCC>`.
- **Cosmetic attributes are not validated**: they do not affect parameterization, but
  they break naive re-reads (see the reload trap) and are not portable beyond `.offxml`.
- **Resource use**: each helper call builds and minimizes an OpenMM `System` on CPU —
  seconds for this ligand. For hundreds of parameter scans, confirm the budget first;
  do not launch large GPU jobs without explicit approval.
- **Never present a modified force field as Sage**: name the file, the base force field
  and version, and every changed parameter in your report.

## Sources

Read at `openforcefield/openff-toolkit` `main` @
`ef074d58a2e7b03625aae7cd3aff87a8ee46c38b` (2026-10-16); byte-identical on `stable` @
`3673715c3223e1d0572a41a346fae113e8aff287` (2026-08-12).

- Notebook: <https://github.com/openforcefield/openff-toolkit/blob/ef074d58a2e7b03625aae7cd3aff87a8ee46c38b/examples/forcefield_modification/forcefield_modification.ipynb>
- Directory: <https://github.com/openforcefield/openff-toolkit/tree/main/examples/forcefield_modification>
- `examples/forcefield_modification/README.md`, `examples/forcefield_modification/.gitignore`,
  `examples/README.md`, `devtools/conda-envs/examples.yaml`
- API signatures verified in `openff/toolkit/typing/engines/smirnoff/parameters.py`
  (`ParameterType.add_cosmetic_attribute`, `ParameterHandler.add_parameter` —
  append-to-end precedence, `allow_duplicate_smirks`, `AngleHandler.AngleType` —
  `angle` in `degree`, `k` in `kilocalorie / mole / degree**2`) and in
  `openff/toolkit/typing/engines/smirnoff/forcefield.py` (`ForceField.__init__`
  `allow_cosmetic_attributes`, `ForceField.to_file` `discard_cosmetic_attributes`,
  `label_molecules`, `create_openmm_system`, `get_partial_charges`)
- Force field values from `openforcefield/openff-forcefields` @
  `be06586b4b8056dac037278a3b717eef65935bce` (`openff-2.3.0.offxml`, ids `a28`/`a28a`)

Upstream licenses: MIT (`openff-toolkit`, `openff-interchange`, `openff-forcefields`).

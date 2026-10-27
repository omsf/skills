---
name: openff-inspect-parameters
description: Inspect which SMIRNOFF parameters a ForceField assigns to atoms, bonds, angles and torsions of specific molecules, and audit parameter coverage over a molecule set, for the openff-toolkit inspect_assigned_parameters example. Use for questions about which SMIRKS matched, parameter IDs (b1/a28/t17), generic vs specific parameter precedence, or get_molecule_parameterIDs coverage reports.
---

# OpenFF Example: Inspecting SMIRNOFF Parameters

## Purpose and scope

Two related tasks with the OpenFF Toolkit:

1. **Per-atom labelling** — for one molecule, list the parameter (with its `id` and
   `smirks`) applied to each bond, angle, torsion, improper torsion and vdW term.
2. **Dataset audit** — for a list of molecules, build the parameter-ID → molecules and
   molecule → parameter-ID maps, to check that a new parameter is actually used or that
   a dataset is fully covered.

Scope ends at *reporting* assignment. This is the diagnostic you run *before* changing
parameters; the editing workflow lives in the `openff-forcefield-modification` skill.

Upstream files: `examples/inspect_assigned_parameters/inspect_assigned_parameters.ipynb`,
`examples/inspect_assigned_parameters/README.md`. A historical, more elaborate coverage
report lives under `examples/deprecated/check_dataset_parameter_coverage/` — treat that
one as legacy (it is excluded from CI) and prefer the two approaches here.

## Required inputs and clarifications

1. **Molecules with explicit hydrogens.** `get_molecule_parameterIDs` documents that it
   wants molecules with explicit hydrogens; the notebook passes explicit-H SMILES. A
   molecule without hydrogens will simply have no H-related parameters to report.
2. **Which force field.** The notebook uses `openff-2.3.0.offxml` (Sage 2.3.0). Parameter
   IDs (`b1`, `a28`, `t17`, ...) are only meaningful relative to a specific force field
   *and version*; always state both when reporting.
3. **Which molecules.** For the audit path, ask for the molecule source (SMILES list,
   SDF, dataset) and whether duplicates are expected — duplicates raise.
4. **Purpose of the report.** "Which parameters does this molecule use?" and "which of
   my new parameters are exercised by this dataset?" need different output shapes; ask
   if unclear.

## Environment setup and version caveats

```bash
mamba create -n openff-toolkit -c conda-forge "openff-toolkit=0.19.0"
mamba activate openff-toolkit
```

- `label_molecules` and `get_molecule_parameterIDs` are pure-Python/SMIRKS matching: no
  OpenMM, no Interchange, no charge model, no NAGL. A minimal `openff-toolkit` install
  (plus RDKit or another toolkit for SMILES parsing) is enough.
- `openff-2.3.0.offxml` resolves from the `openforcefields` package, not from
  `openff/toolkit/data/`.
- The notebook imports `Topology` but only uses it in Example 1; `label_molecules` takes
  a `Topology`, not a bare `Molecule`.
- Aromaticity: matching depends on the force field's `aromaticity_model` (Sage uses
  `OEAroModel_MDL`) and the toolkit's perception of the input molecule. If you build
  molecules with a non-default aromaticity model you will see
  `IncompatibleAromaticityError` / `AromaticityModelError` rather than wrong output.
- Verified against `openff-toolkit` 0.19.0, RDKit 2026.9.1, `openff-2.3.0.offxml`
  sha256 `7a0d7195a4b717e29fe14aa2cf1bc3203e87ce4aef9cc2ef10fb5145b7edcf63`.

## Steps

### 1. Label one molecule

```python
from openff.toolkit import ForceField, Molecule, Topology

molecule = Molecule.from_smiles("C=C(N)(C)")
topology = Topology.from_molecules([molecule])

forcefield = ForceField("openff-2.3.0.offxml")
molecule_force_list = forcefield.label_molecules(topology)

for mol_idx, mol_forces in enumerate(molecule_force_list):
    print(f"Forces for molecule {mol_idx}")
    for force_tag, force_dict in mol_forces.items():
        print(f"\n{force_tag}:")
        for atom_indices, parameter in force_dict.items():
            atomstr = "".join(f"{idx:>3}" for idx in atom_indices)
            print(f"atoms: {atomstr}  parameter_id: {parameter.id}  smirks {parameter.smirks}")
```

`label_molecules` returns **one dict per molecule in the topology**; each dict maps a
handler tag (`Constraints`, `Bonds`, `Angles`, `ProperTorsions`, `ImproperTorsions`,
`vdW`, `Electrostatics`, `LibraryCharges`, `NAGLCharges`) to a dict keyed by tuples of
**0-indexed atom indices within that molecule** with `ParameterType` values. Keys are
`ValenceDict`/`ImproperDict` for the valence handlers, so atom order inside a key is
normalized for torsions and impropers.

`Electrostatics`, `LibraryCharges` and `NAGLCharges` can come back as **empty
dicts** (or be absent, depending on the force field). `Constraints` does have
SMIRKS matches in this example; inspect it like the other parameterized terms. Do not read an empty dict as "this term is missing".

### 2. Work with the labels programmatically

The notebook derives atom indices rather than hardcoding them, because SMILES atom
ordering is not guaranteed to match your reference:

```python
labels = forcefield.label_molecules(molecule.to_topology())[0]

# Which parameters apply to a specific bond?
for atoms, parameter in labels["Bonds"].items():
    if parameter.id == "b6":
        print(atoms, parameter.smirks)

# Group by parameter to see how often each is used
from collections import Counter

usage = {tag: Counter(p.id for p in d.values()) for tag, d in labels.items()}
print(usage["Angles"])
```

### 3. Audit a molecule set

```python
from openff.toolkit import ForceField, Molecule
from openff.toolkit.utils import get_molecule_parameterIDs

molecule_smiles = [
    "[H][C@@]1([C@](OC1([H])[H])([H])O[H])O[H]",
    "[H]C1(C(C(C1([H])O[H])([H])[H])([H])O[H])[H]",
    "[H]C1(C(C(O1)([H])[H])([H])O[H])[H]",
    "[H][C@@]1(C([C@@](OC1([H])[H])([H])O[H])([H])[H])O[H]",
    "[H][C@@]1(C(O[C@](O1)([H])O[H])([H])[H])O[H]",
    "[H][C@]1(C([C@](O1)([H])O[H])([H])[H])O[H]",
    "[H][C@]1(C(OC(O1)([H])[H])([H])[H])O[H]",
    "[H][C@@]1([C@](C(C(C1([H])[H])([H])[H])([H])[H])([H])O[H])O[H]",
    "[H][C@]1([C@@](C(OC(C(C1([H])[H])([H])[H])([H])[H])([H])[H])([H])O[H])O[H]",
    "[H][C@@]1([C@](OC(C(C1([H])[H])([H])[H])([H])[H])([H])O[H])O[H]",
]

mols = [Molecule.from_smiles(smiles) for smiles in molecule_smiles]
forcefield = ForceField("openff-2.3.0.offxml")

parameters_by_molecule, parameters_by_ID = get_molecule_parameterIDs(mols, forcefield)

print("Molecules with parameter IDs:")
for pid in sorted(parameters_by_ID):
    print(pid)
    for smiles in parameters_by_ID[pid]:
        print("\t", smiles)
```

Return contract (from the docstring, confirmed by running it):

- `parameters_by_molecule`: `{isomeric_smiles: [parameter_id, ...]}` — the list is **not
  deduplicated**; a parameter used 5 times appears 5 times.
- `parameters_by_ID`: `{parameter_id: {isomeric_smiles, ...}}` — a `set`, so no
  frequency information.

### 4. Turn the audit into a coverage check

```python
# New parameters you just added, by id
new_parameter_ids = {"a99"}          # placeholder: substitute your real ids
covered = new_parameter_ids & set(parameters_by_ID)
print("new parameters exercised:", covered or "none")

uncovered = new_parameter_ids - set(parameters_by_ID)
if uncovered:
    print("WARNING: these new parameters never matched:", uncovered)
```

This is the pattern the notebook is advertising: adding a parameter and verifying it is
actually assigned to the dataset it was fitted for.

## Validation and expected outputs

- Example 1 for `C=C(N)(C)` with Sage 2.3.0 prints 9 handlers. Observed: 7 `Constraints`
  entries all `id: c1 / [#1:1]-[*:2]`; 10 `Bonds` (e.g. `0 1 -> b6
  [#6X3:1]=[#6X3:2]`, `1 2 -> b8 [#6X3:1]-[#7X3:2]`, `0 4 -> b85
  [#6X3:1]-[#1:2]`); 15 `Angles` (e.g. `4 0 5 -> a12 [#1:1]-[#6X3:2]-[#1:3]`,
  `8 3 9 -> a2 [#1:1]-[#6X4:2]-[#1:3]`); 17 `ProperTorsions`; 2 `ImproperTorsions`
  (`i1`, `i4`); and empty `Electrostatics`, `LibraryCharges`, `NAGLCharges` dicts.
  Counts depend on the molecule; the *shape* and the presence of empty dicts do not.
- Example 2 on the ten explicit-H sugars yields **33 distinct parameter IDs** across
  `parameters_by_ID` with Sage 2.3.0.
- Pass criteria for an audit: every molecule appears in `parameters_by_molecule`, every
  parameter you care about appears in `parameters_by_ID`, and no parameter is reported
  for a molecule you did not supply.

## Failure modes and safety

- **`label_molecules` does not fail on missing parameters (verified).** For
  `C[SiH3]` with an unconstrained Sage force field it happily returned a `Bonds` dict
  containing only the three C-H bonds and silently omitted the C-Si and Si-H bonds.
  Only `create_openmm_system`/`create_interchange` raise
  (`UnassignedBondError`, `UnassignedAngleParameterException`,
  `UnassignedValenceParameterException`). **Always compare the number of labelled terms
  against the molecule's actual bond/angle/torsion counts** if you are using this to
  judge coverage:

```python
labels = forcefield.label_molecules(molecule.to_topology())[0]
assert len(labels["Bonds"]) == molecule.n_bonds, "some bonds were not parameterized"
```

- **Duplicates raise.** `get_molecule_parameterIDs` compares isomeric SMILES and raises
  `ValueError("...contains some duplicates: {...}")` if two inputs are the same species.
  Deduplicate first.
- **Keys are isomeric SMILES, not names.** If you need to map results back to your own
  identifiers, keep a SMILES → your-id lookup; `Molecule.name` is not used.
- **Atom indices are per-molecule and 0-indexed**, and the dict keys for torsions are
  order-normalized. Do not assume the tuple order matches how you enumerated bonds.
- **Parameter IDs are force-field-version-specific.** A report that says "`a28`" is
  meaningless without "Sage 2.3.0".
- **`get_molecule_parameterIDs` is O(molecules) SMIRKS matching.** It is fine for
  thousands of molecules but is single-threaded; for very large sets, tell the user the
  expected runtime rather than starting a multi-hour run unannounced.

## Sources

Read at `openforcefield/openff-toolkit` `main` @
`ef074d58a2e7b03625aae7cd3aff87a8ee46c38b` (2026-10-16); byte-identical on `stable` @
`3673715c3223e1d0572a41a346fae113e8aff287` (2026-08-12).

- Notebook: <https://github.com/openforcefield/openff-toolkit/blob/ef074d58a2e7b03625aae7cd3aff87a8ee46c38b/examples/inspect_assigned_parameters/inspect_assigned_parameters.ipynb>
- Directory: <https://github.com/openforcefield/openff-toolkit/tree/main/examples/inspect_assigned_parameters>
- `examples/inspect_assigned_parameters/README.md`, `examples/README.md`,
  `devtools/conda-envs/examples.yaml`, `.github/workflows/examples.yml`
- Historical alternative (excluded from CI): `examples/deprecated/check_dataset_parameter_coverage/check_parameter_coverage.ipynb`
  <https://github.com/openforcefield/openff-toolkit/tree/main/examples/deprecated/check_dataset_parameter_coverage>
- `ForceField.label_molecules` (return contract) in
  `openff/toolkit/typing/engines/smirnoff/forcefield.py`;
  `get_molecule_parameterIDs` in `openff/toolkit/utils/utils.py` (line ~704)
- Force field `openff-2.3.0.offxml` from `openforcefield/openff-forcefields` @
  `be06586b4b8056dac037278a3b717eef65935bce`

Upstream licenses: MIT (`openff-toolkit`, `openff-forcefields`).

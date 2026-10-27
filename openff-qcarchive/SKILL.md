---
name: openff-qcarchive
description: Pull records out of QCArchive with the qcportal client and turn them into OpenFF Molecule objects (and back) via Molecule.from_qcschema / to_qcschema, including CMILES-based bond-order recovery, atom-order verification and a QCEngine energy — for the openff-toolkit QCArchive_interface example. Use for requests about QCArchive, QCFractal, MolSSI datasets, QCSchema, or reusing public QM geometries in OpenFF workflows.
---

# OpenFF Example: QCArchive Interface

## Purpose and scope

Connect to a QCFractal/QCArchive instance, fetch an optimization dataset entry, and
convert it into an OpenFF `Molecule` with correct bond orders, formal charges and
stereochemistry by way of the CMILES field
`canonical_isomeric_explicit_hydrogen_mapped_smiles`. Then verify the atom ordering
against the raw QCArchive geometry, optionally run a QCEngine calculation, and convert
the molecule back to QCSchema.

Scope ends at getting a trustworthy `Molecule` (plus one geometry) out of QCArchive.
Electronic structure is not computed here beyond a toy RDKit UFF single point.

Upstream file: `examples/QCArchive_interface/QCarchive_interface.ipynb`. There is no
README in that directory.

## Required inputs and clarifications

1. **Which instance.** The example uses the public instance
   `https://api.qcarchive.molssi.org:443`. Any other host, credentials or
   `PortalClient` options (e.g. a Fractal server URL) must come from the user.
2. **Which dataset and entry.** The example uses
   `dataset_type="optimization"`, `dataset_name="Kinase Inhibitors: WBO Distributions"`
   and `dataset.entry_names[-1]`. Dataset names are not stable API — confirm the name
   against `client.list_datasets_table()` before assuming it exists.
3. **Which geometry you want.** For an optimization entry the example uses
   `entry.initial_molecule`; the *optimized* geometry is a separate record. Ask.
4. **Network and data volume — get approval.** Every step here is a network call to a
   public server. `list_datasets_table()` is cheap; `get_dataset(...)` on a large set
   (e.g. `OpenFF VEHICLe Set 1`, 25500 records) is not. Do not download bulk records
   without explicit user approval. `qcengine` is only needed for step 5 — skip it if the
   user just wants `Molecule` objects.

## Environment setup and version caveats

```bash
mamba create -n openff-toolkit -c conda-forge "openff-toolkit=0.19.0" qcportal qcelemental qcengine networkx
mamba activate openff-toolkit
```

- The examples environment `devtools/conda-envs/examples.yaml` on main uses `qcportal >=0.70` (stable uses
  `>=0.50`),
  `qcelemental`, `qcengine`, `pydantic =2`, `networkx`, `rdkit`.
- Verified against `openff-toolkit` 0.19.0, `qcportal` 0.71, `qcelemental` 0.51.4,
  `qcengine` 0.52.0, RDKit 2026.9.1.
- **Two versions of this notebook exist and they are not interchangeable.** On branch
  `stable` it uses `entry.dict()` and a QCSchema **v1** QCEngine task; on `main` it uses
  `entry.model_dump()` and a QCSchema **v2** task (`schema_name:
  "qcschema_atomic_input"`, `specification: {...}`).
  - `Molecule.to_qcschema()` returns `qcelemental.models.v2.Molecule` only on `main`
    (PR #2238, unreleased in 0.19.0). On the **released** toolkit it returns a **v1**
    molecule, so the `main` notebook's v2 task fails with
    `QCEngine Input Error: Error creating 'AtomicInput' ... molecule Input should be a
    valid dictionary or instance of Molecule` — reproduced with `qcengine` 0.52.0 /
    `qcelemental` 0.51.4. Use the v1 task form below unless you have a toolkit build
    whose `to_qcschema()` emits v2.
  - `entry.dict()` works with Pydantic 2 but emits
    `PydanticDeprecatedSince20` (removed in Pydantic 3). Prefer `entry.model_dump()`
    for the dict argument; it is equivalent for `from_qcschema`.
- `client.get_molecules(...)` returns a `qcelemental` Molecule (in this stack a
  `qcelemental.models._v1v2.molecule.Molecule`), **not** an OpenFF `Molecule`.

## Steps

### 1. Connect and pick a record

```python
import qcportal

client = qcportal.PortalClient("https://api.qcarchive.molssi.org:443")

# Cheap: prints a table of dataset id, type, record_count, name
print(client.list_datasets_table())

dataset = client.get_dataset(
    dataset_type="optimization",
    dataset_name="Kinase Inhibitors: WBO Distributions",
)
entry = dataset.get_entry(entry_name=dataset.entry_names[-1])
print(entry.name)
print(sorted(entry.model_dump().keys()))   # ['additional_keywords', 'attributes',
                                           #  'comment', 'initial_molecule', 'name']
```

### 2. Build a `Molecule` from the entry

```python
from openff.toolkit import Molecule

molecule_from_entry = Molecule.from_qcschema(entry)
molecule_from_dict = Molecule.from_qcschema(entry.model_dump())
assert molecule_from_entry == molecule_from_dict
molecule = molecule_from_entry
```

`from_qcschema` accepts an entry, a QCSchema molecule, or a dict of either. It recovers
the Kekulé structure from `canonical_isomeric_explicit_hydrogen_mapped_smiles`, looked
up in this order: the entry's `attributes`, then `identifiers`, then `extras`. A raw
QCSchema molecule with no CMILES field cannot be converted — that is exactly why the
method exists.

### 3. Verify atom ordering against the raw record

This is the step that makes the conversion trustworthy; do not skip it.

```python
initial_molecule = client.get_molecules(entry.initial_molecule.id)
print(type(initial_molecule))   # qcelemental Molecule, not OpenFF

assert molecule.n_atoms == len(initial_molecule.atomic_numbers), "atom count mismatch"
for off_atom, atomic_number in zip(molecule.atoms, initial_molecule.atomic_numbers):
    assert off_atom.atomic_number == atomic_number, "atom order mismatch"

print(molecule.to_hill_formula(), initial_molecule.get_molecular_formula())
```

```python
import networkx as nx

initial_network = nx.Graph()
for index, atomic_number in enumerate(initial_molecule.atomic_numbers):
    initial_network.add_node(index, atomic_number=atomic_number)
for bond in initial_molecule.connectivity:
    initial_network.add_edge(*bond[:2])

isomorphic, atom_map = Molecule.are_isomorphic(
    molecule,
    initial_network,
    return_atom_map=True,
    aromatic_matching=False,
    formal_charge_matching=False,
    bond_order_matching=False,
    bond_stereochemistry_matching=False,
    atom_stereochemistry_matching=False,
)
assert isomorphic
assert all(i == j for i, j in atom_map.items()), "graph is isomorphic but permuted"
```

The `*_matching=False` flags are deliberate: the QCArchive graph carries only
connectivity, so any stricter comparison would fail. `are_isomorphic` accepts either a
`Molecule` or a `networkx.Graph`.

### 4. Convert back to QCSchema

```python
qc_molecule = molecule.to_qcschema()
print(qc_molecule.schema_name, qc_molecule.schema_version)   # qcschema_molecule 2
print(qc_molecule.extras["canonical_isomeric_explicit_hydrogen_mapped_smiles"])
```

`to_qcschema` requires at least one conformer and raises `InvalidConformerError`
otherwise; pass `conformer=<index>` to select one. It writes the mapped SMILES into
both `extras` and `identifiers`.

### 5. Optional: a QCEngine single point

```python
import qcengine

rdkit_task = {
    "schema_name": "qcschema_input",          # QCSchema v1 task form
    "molecule": qc_molecule,
    "driver": "energy",
    "model": {"method": "uff", "basis": None},
    "keywords": {"scf_type": "df"},
}
result = qcengine.compute(rdkit_task, "rdkit")
assert result.success, result.error
print(result.return_result)                   # hartree
```

Always check `result.success` (or `isinstance(result, FailedOperation)`) before reading
`return_result`; a failed computation returns a `FailedOperation`, not a result object.

## Validation and expected outputs

- `client.list_datasets_table()` returns a **string** (a rendered table), not a
  dataframe — index/print it, do not call `.head()`.
- Reproduced on the live public instance: `Kinase Inhibitors: WBO Distributions` is an
  `optimization` dataset; the last entry is
  `cs(=o)(=o)ccncc1ccc(o1)c2ccc3c(c2)c(ncn3)nc4ccc(c(c4)cl)occ5cccc(c5)f-99`
  with 66 atoms, 1 conformer, formula `C29H26ClFN4O4S`.
- `molecule.properties` after `from_qcschema` on an optimization entry is
  `{"initial_molecules": {0: <qc molecule id>}}` — conformer index → QCArchive id.
- The atom-order assertion and the identity `atom_map` assertion must both pass. A
  non-identity `atom_map` means the CMILES mapping and the stored geometry disagree;
  stop and investigate rather than using the geometry.
- `qcengine.compute(rdkit_task, "rdkit")` on that entry returned
  `return_result = 0.05930479138457709` hartree with `success=True` (UFF, RDKit
  backend). The absolute number is a toy UFF energy; treat it as a plumbing check.
- A round trip through `Molecule.to_qcschema()` → `Molecule.from_qcschema(...)` must
  reproduce the same molecule (`Molecule.from_qcschema(molecule.to_qcschema()) == molecule`).

## Failure modes and safety

- **Missing CMILES**: `from_qcschema` on a molecule record with no
  `canonical_isomeric_explicit_hydrogen_mapped_smiles` in `attributes`/`identifiers`/
  `extras` raises; QCArchive geometries alone do not encode bond orders or formal
  charges.
- **Undefined stereochemistry**: `from_qcschema` raises unless you pass
  `allow_undefined_stereo=True`. Do not flip that flag to make an error go away — it
  changes the chemistry you get back.
- **Version mismatch (the big one)**: on a released toolkit, the `main`-branch QCEngine
  task form fails at input validation. If you hit `Error creating 'AtomicInput'`, check
  `molecule.to_qcschema().schema_version` and use the matching task schema.
- **Dataset drift**: dataset names, entry names and record counts on the public instance
  change over time. Never hardcode a record count or assume `entry_names[-1]` is the
  same molecule it was upstream.
- **Shared public server**: the MolSSI instance is a shared resource. Cache downloaded
  records to disk, avoid re-querying in loops, and prefer `client.query_molecules(...)`
  with a `limit` for spot checks.
- **`qcengine` may shell out**: some QC backends run external binaries. Confirm which
  backend is intended; here it is RDKit/UFF, which is in-process and cheap.

## Sources

Read at `openforcefield/openff-toolkit` `stable` @
`3673715c3223e1d0572a41a346fae113e8aff287` (2026-08-12) and `main` @
`ef074d58a2e7b03625aae7cd3aff87a8ee46c38b` (2026-10-16); the notebook differs between
the two branches (v1 vs v2 QCSchema task, `dict()` vs `model_dump()`).

- Notebook (stable, v1 task): <https://github.com/openforcefield/openff-toolkit/blob/3673715c3223e1d0572a41a346fae113e8aff287/examples/QCArchive_interface/QCarchive_interface.ipynb>
- Notebook (main, v2 task): <https://github.com/openforcefield/openff-toolkit/blob/ef074d58a2e7b03625aae7cd3aff87a8ee46c38b/examples/QCArchive_interface/QCarchive_interface.ipynb>
- Directory: <https://github.com/openforcefield/openff-toolkit/tree/stable/examples/QCArchive_interface>
- `Molecule.from_qcschema` / `Molecule.to_qcschema` docstrings and CMILES lookup order:
  `openff/toolkit/topology/molecule.py` (main: `to_qcschema` returns
  `qcel.models.v2.Molecule`; stable: `qcel.models.Molecule`). Release note for the v2
  change: `docs/releasehistory.md` (PR #2238)
- `devtools/conda-envs/examples.yaml`, `.github/workflows/examples.yml`
- QCFractal dataset user guide: <https://docs.qcarchive.molssi.org/user_guide/datasets/index.html>;
  QCEngine: <https://molssi.github.io/QCEngine/>

Upstream licenses: MIT (`openff-toolkit`); `qcportal`/`qcelemental`/`qcengine` are
BSD-3-Clause (MolSSI).

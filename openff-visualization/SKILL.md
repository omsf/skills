---
name: openff-visualization
description: Render OpenFF Molecule and Topology objects in Jupyter with the toolkit's built-in visualization backends (RDKit, OpenEye, NGLview), including implicit display, explicit Molecule.visualize()/Topology.visualize() calls, hydrogen and sizing options, and NGLview troubleshooting. Use for requests about visualizing molecules/topologies with OpenFF, or the openff-toolkit visualization example.
---

# Visualizations with OpenFF

## Purpose and scope

Show how `openff.toolkit.topology.Molecule` and `Topology` render in Jupyter — both
implicitly (evaluating the object in a cell) and explicitly via `.visualize()` — and
how to choose and configure a backend. Upstream files:
`examples/visualization/visualization.ipynb` and `examples/visualization/README.md`.

Scope is display only. This is not a structure-preparation guide; the toolkit's
rendering is a 2D depiction (RDKit/OpenEye) or a 3D widget (NGLview), not a
publication-quality figure pipeline.

## Required inputs and clarifications

1. **Display target**: a `Molecule` (SMILES/file/`Molecule`s from a `Topology`) or a
   whole `Topology` (e.g. a PDB). Ask which, since the two APIs differ.
2. **3D or 2D**: NGLview requires conformers/positions; RDKit and OpenEye are 2D
   depictions and need no coordinates but cannot show stereochemistry-in-3D.
3. **OpenEye availability**: the OpenEye backend needs an installed *and licensed*
   OpenEye toolkit. Do not assume it is present.
4. **Notebook vs headless/script**: these APIs return widget/`SVG`/`Image` objects
   intended for IPython display. Confirm whether the user is actually in a Jupyter
   environment before recommending NGLview troubleshooting.

## Environment setup and version caveats

```bash
mamba create -n openff-toolkit -c conda-forge "openff-toolkit=0.19.0"
mamba activate openff-toolkit
# 3D widget backend (optional):
mamba install -c conda-forge nglview
```

- RDKit is installed automatically with the `openff-toolkit` package (not with
  `openff-toolkit-base`) and provides the 2D default backend.
- The examples environment (`devtools/conda-envs/examples.yaml`, name
  `openff-toolkit-examples`) adds `nglview` and `mdtraj`; CI runs the notebook with
  `pytest --nbval-lax --ignore=examples/deprecated examples`.
- NGLview is the fragile dependency. The older instructions
  `jupyter-nbextension enable nglview --py --sys-prefix` and
  `jupyter labextension install nglview-js-widgets` /
  `jupyter-labextension install @jupyter-widgets/jupyterlab-manager` **predate
  JupyterLab 4**; `jupyter labextension install` no longer exists there. With
  NGLview ≥ 3.0.0 and a modern Jupyter stack you normally need no extension step, but
  verify the widget actually renders before promising it.
- OpenEye installation (`mamba install -c openeye openeye-toolkits`) also requires a
  license file; `oechem.OEChemIsLicensed()` is the guard the notebook uses.

## Steps

### 1. Implicit visualization

Evaluating a `Molecule` as the last expression in a cell displays it. The toolkit's
`_ipython_display_` tries backends in order: **NGLview → RDKit → OpenEye**. Because
molecules without conformers cannot use NGLview, a SMILES-only molecule falls through
to the RDKit 2D depiction; once conformers exist, NGLview takes over.

```python
from IPython.display import display
from openff.toolkit import Molecule

m = Molecule.from_smiles("CCCCOCC")
m            # implicit: NGLview if conformers exist, else RDKit, else OpenEye
display(m)   # equivalent explicit display call
```

`Topology` has no implicit representation (only `__repr__`/`_repr_pretty_`); you must
call `Topology.visualize()` explicitly.

### 2. Explicit `Molecule.visualize()`

Signature: `visualize(backend="rdkit", width=500, height=300, show_all_hydrogens=True)`.

```python
m.visualize()                                   # default backend is "rdkit"
m.visualize(backend="rdkit")
m.visualize(backend="rdkit", show_all_hydrogens=False)
m.visualize(width=1000, height=100)             # ignored with backend="nglview"

m.generate_conformers()
m.visualize(backend="nglview")                  # 3D widget with a trajectory player
```

- `width`/`height` apply only to the RDKit/OpenEye backends. Passing non-default
  values with `backend="nglview"` emits a `UserWarning`.
- `show_all_hydrogens=False` maps to `Chem.RemoveHs(...)` for RDKit and to an NGL
  selection of `"NOT hydrogen"`.
- A `Molecule` with several conformers gets a trajectory player in NGLview; a rigid
  molecule (e.g. benzene) will not show one.

OpenEye backend, guarded:

```python
try:
    from openeye import oechem

    assert oechem.OEChemIsLicensed()
    m.visualize(backend="openeye")
except (ImportError, AssertionError):
    print('Visualizing with `backend="openeye"` requires the OpenEye Toolkits')
```

### 3. Explicit `Topology.visualize()`

```python
from openff.toolkit import Topology
from openff.toolkit.utils import get_data_file_path

topology = Topology.from_pdb(
    get_data_file_path("systems/test_systems/T4_lysozyme_water_ions.pdb")
)
w = topology.visualize()
w   # display the NGLWidget
```

`Topology.visualize()` takes no `backend` argument and always uses NGLview. It requires
positions (`MissingConformersError` otherwise) and adds representations for water
(`line`), ions (`spacefill`), protein (`cartoon`), and everything else (`licorice`),
plus a unit-cell representation.

## Validation and expected outputs

- RDKit backend returns an `IPython.display.SVG`; OpenEye returns an
  `IPython.display.Image`; NGLview returns an `nglview.NGLWidget` that renders
  interactively in a Jupyter cell.
- Quick non-visual smoke check in a headless script:

```python
from IPython.display import SVG

m.visualize()                       # returns an SVG object even outside a notebook
assert isinstance(m.visualize(), SVG)
```

- For a `Topology`, assert positions exist before visualizing:
  `assert topology.get_positions() is not None`, then check the returned widget's
  particle count equals the topology's particle count.

## Failure modes and safety

- **No conformers with NGLview**: `.visualize(backend="nglview")` raises
  `MissingConformersError`, which subclasses `ValueError`; call
  `Molecule.generate_conformers()` first (an ETKDG/OpenEye conformer generation step,
  so it costs time and is stochastic).
- **`ensure_correct_connectivity=True` is not implemented**:
  `Topology.visualize(ensure_correct_connectivity=True)` raises
  `ValueError("`ensure_correct_connectivity` not (yet) implemented ...")`. The default
  (positions-inferred) connectivity is documented as possibly inaccurate, and this is a
  deliberate performance trade-off — do not present it as a bug.
- **Silent widget failure**: NGLview can import successfully yet render nothing if the
  Jupyter front-end extension/labextension is missing. That is a front-end problem, not
  a toolkit error; do not start editing chemistry code to fix it.
- **OpenEye confusion**: `backend="openeye"` on a machine without a license raises
  inside OpenEye; guard with `OEChemIsLicensed()` as above rather than catching broadly.
- **Do not compare depictions for scientific claims**: 2D RDKit/OpenEye output does not
  validate stereochemistry or conformers. Check the underlying `Molecule` properties.
- **Large topologies**: `Topology` visualization is explicitly not automatic because
  generating NGLview structures for large systems is expensive. Ask before rendering
  solvated boxes or big complexes, and prefer selection-specific views.

## Sources

Read at commit `ef074d58a2e7b03625aae7cd3aff87a8ee46c38b`
(`openforcefield/openff-toolkit`, `main`, 2026-10-16); `examples/visualization/` is
byte-identical on branch `stable` at `3673715c3223e1d0572a41a346fae113e8aff287`.

- Notebook: `examples/visualization/visualization.ipynb`
  <https://github.com/openforcefield/openff-toolkit/blob/ef074d58a2e7b03625aae7cd3aff87a8ee46c38b/examples/visualization/visualization.ipynb>
- Directory: <https://github.com/openforcefield/openff-toolkit/tree/stable/examples/visualization>
- `examples/visualization/README.md` (backend list and the NGLview setup commands)
- API: `Molecule.visualize` / `_ipython_display_` in
  `openff/toolkit/topology/molecule.py`; `Topology.visualize` in
  `openff/toolkit/topology/topology.py`; `MissingConformersError` and
  `MissingOptionalDependencyError` in `openff/toolkit/utils/exceptions.py`
- Docs: `docs/installation.md` (RDKit/AmberTools by default, OpenEye optional and licensed)
- Data: `openff/toolkit/data/systems/test_systems/T4_lysozyme_water_ions.pdb`

Licenses: MIT (`openforcefield/openff-toolkit`). NGLview, RDKit, and OpenEye carry
their own licenses; OpenEye requires a (free for some uses) license file.

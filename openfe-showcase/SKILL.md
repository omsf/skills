---
name: openfe-showcase
description: Set up and analyse an OpenFE relative binding free energy calculation in Python using the TYK2 showcase (components, Kartograf/LOMAP atom mappings, ligand networks, ChemicalSystems, RelativeHybridTopologyProtocol settings, Transformation create/dry-run/to_json, and gather-based analysis). Use when the CLI is too rigid, when a user wants to understand or customise OpenFE protocol settings in Python, when working through the openfe_showcase.ipynb notebook, or when they need the mapping/network objects that the CLI hides.
---

# OpenFE Showcase (Python API)

## Purpose and scope

The showcase notebook walks the full TYK2 RBFE setup in Python: load ligands,
build atom mappings with two mappers, build and visualise three network
topologies, define `ChemicalSystem`s, inspect and override protocol settings,
create `Transformation`s, dry-run the OpenMM hybrid-topology protocol, write
transformations to JSON for `openfe quickrun`, and gather free energies.

Use this when the CLI (`openfe plan-rbfe-network`) cannot express what you need —
custom settings, single-edge calculations, or direct access to mapping/network
objects. For the plain CLI path see `openfe-rbfe-cli`. Scope: RBFE only; OpenFE
also ships ABFE, AHFE, SepTop, and MD protocols, which the notebook links to but
does not demo.

## Required inputs and clarifications

- **System**: two ligands plus a protein. Tutorial files are
  `inputs/tyk2_ligands.sdf` (10 ligands) and `inputs/tyk2_protein.pdb`
  (capped and protonated TYK2).
- **Which ligand pair?** The notebook hard-codes `lig_ejm_31 -> lig_ejm_47`.
- **Protocol settings**: which equilibration/production lengths and platform?
  See the cost warning in step 5.
- **Approval.** Even the notebook's shortened production run is GPU work, and the
  Colab bootstrap installs a large conda environment. Get user approval before
  executing any simulation cell.

## Environment

Notebook pins (ExampleNotebooks `environment.yaml`):

```bash
micromamba create -n openfe-notebooks -c conda-forge \
    python=3.13 'openfe==1.13.*' 'cinnabar>=0.6.0' 'kartograf~=2.0' \
    'konnektor~=0.5.0' 'openfe-analysis>=0.5.0' 'lomap2>=3.2.1' \
    'nglview<4' 'setuptools<80' numpy notebook ipywidgets
micromamba activate openfe-notebooks
```

That file *also* pip-installs `git+https://github.com/OpenFreeEnergy/openfe@main`,
which supersedes the `1.13.*` conda pin — notebook CI tracks `main`. Drop the pip
line to avoid a moving development install. This is still not a complete lockfile;
record resolved versions. APIs below were source-checked on main, not run-checked
against the release. Confirm installed APIs before adopting main-only features.

Optional extras: `py3Dmol` (needed by `LigandAtomMapping.view_3d()`),
`matplotlib` + `networkx` + `rdkit` (needed by `plot_atommapping_network`).

The notebook's Colab cell installs the **OpenFEforge 1.7.0** single-file
installer via `condacolab`. That is historical and pins OpenFE 1.7.0 — do not
reuse it for new work; skip it locally and follow the current install docs.

Get the inputs (identical to the RBFE tutorial files):

```bash
mkdir -p showcase_work/inputs && cd showcase_work/inputs
openfe fetch rbfe-tutorial
ls   # tyk2_ligands.sdf  tyk2_protein.pdb
cd ..   # all Python paths below are relative to showcase_work, not inputs/
```

## Procedure

### 1. Ligands as `SmallMoleculeComponent`s

```python
from rdkit import Chem
from openfe import SmallMoleculeComponent

ligand_mols = [SmallMoleculeComponent(m)
               for m in Chem.SDMolSupplier("inputs/tyk2_ligands.sdf", removeHs=False)]

# equivalent route via the OpenFF toolkit
from openff.toolkit import Molecule
ligand_mols = [SmallMoleculeComponent.from_openff(m)
               for m in Molecule.from_file("inputs/tyk2_ligands.sdf")]

print(ligand_mols[0].name)                   # lig_ejm_31  (from the SDF _Name)
print(ligand_mols[0].to_openff().__class__)  # openff.toolkit.topology.Molecule
```

### 2. Atom mappings

```python
from openfe.setup import KartografAtomMapper, LomapAtomMapper

mapping_k = next(KartografAtomMapper(atom_map_hydrogens=True)
                 .suggest_mappings(ligand_mols[0], ligand_mols[4]))
mapping_l = next(LomapAtomMapper()
                 .suggest_mappings(ligand_mols[0], ligand_mols[4]))

mapping_k                       # 2D representation (Jupyter repr)
mapping_k.view_3d(show_atomIDs=True)   # needs py3Dmol
mapping_k.draw_to_file("tyk2_edge.png")  # RDKit PNG
```

`suggest_mappings` is a generator; `next()` raises `StopIteration` if no mapping
exists. `view_3d` and `draw_to_file` live on `gufe.LigandAtomMapping`. In 3D,
shared sphere colours are mapped atoms; uncoloured atoms become dummy atoms in
the opposite end state.

### 3. Ligand networks

```python
import openfe
from openfe.setup import KartografAtomMapper
from openfe.setup.ligand_network_planning import (
    generate_minimal_spanning_network, generate_lomap_network, generate_radial_network)
from openfe.utils.atommapping_network_plotting import plot_atommapping_network

mappers = [KartografAtomMapper()]
scorer = openfe.lomap_scorers.default_lomap_score

mst_network = generate_minimal_spanning_network(ligands=ligand_mols, scorer=scorer, mappers=mappers)
lomap_network = generate_lomap_network(ligands=ligand_mols, scorer=scorer, mappers=mappers)
radial_network = generate_radial_network(
    ligands=ligand_mols[1:], central_ligand=ligand_mols[0], mappers=mappers)

fig = plot_atommapping_network(mst_network)   # returns a matplotlib Figure

edge = list(mst_network.edges)[1]
print(edge.componentA.smiles, edge.componentB.smiles)
print(edge.componentA_to_componentB)          # {indexA: indexB}

with open("network_store.graphml", "w") as fh:
    fh.write(mst_network.to_graphml())
```

`generate_lomap_network` is not implemented in OpenFE — it is re-exported from
the external `lomap` package (`lomap2>=3.2.1` in the notebook env). Transformations
are always directed `componentA -> componentB`.

### 4. Components and `ChemicalSystem`s

```python
from openfe import ChemicalSystem, ProteinComponent, SolventComponent
from openff.units import unit

protein = ProteinComponent.from_pdb_file("inputs/tyk2_protein.pdb")
solvent = SolventComponent(positive_ion="Na", negative_ion="Cl",
                           neutralize=True, ion_concentration=0.15 * unit.molar)

matches = [e for e in mst_network.edges
           if e.componentA.name == "lig_ejm_31" and e.componentB.name == "lig_ejm_47"]
if len(matches) != 1:
    raise ValueError("Expected directed edge 31 -> 47 is absent; inspect the network")
edge = matches[0]

a_complex = ChemicalSystem({"ligand": edge.componentA, "solvent": solvent, "protein": protein},
                           name=edge.componentA.name)
a_solvent = ChemicalSystem({"ligand": edge.componentA, "solvent": solvent},
                           name=edge.componentA.name)
b_complex = ChemicalSystem({"ligand": edge.componentB, "solvent": solvent, "protein": protein},
                           name=edge.componentB.name)
b_solvent = ChemicalSystem({"ligand": edge.componentB, "solvent": solvent},
                           name=edge.componentB.name)
```

An RBFE edge needs all four systems. `SolventComponent` is a *recipe* — water and
ions are added later by the protocol's solvation settings, not here.

### 5. Protocol and settings

`RelativeHybridTopologyProtocol.default_settings()` returns eleven settings
groups: `protocol_repeats`, `simulation_settings`, `output_settings`,
`alchemical_settings`, `engine_settings`, `integrator_settings`,
`lambda_settings`, `forcefield_settings`, `thermo_settings`,
`solvation_settings`, `partial_charge_settings`.

```python
from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol
from openff.units import unit

# SMOKE TEST ONLY -- do not use these lengths for science.
solvent_settings = RelativeHybridTopologyProtocol.default_settings()
solvent_settings.simulation_settings.equilibration_length = 10 * unit.picosecond
solvent_settings.simulation_settings.production_length = 50 * unit.picosecond
solvent_settings.engine_settings.compute_platform = None   # auto-select CUDA/OpenCL/CPU

complex_settings = RelativeHybridTopologyProtocol.default_settings()
complex_settings.simulation_settings.equilibration_length = 10 * unit.picosecond
complex_settings.simulation_settings.production_length = 50 * unit.picosecond
complex_settings.solvation_settings.solvent_padding = 1 * unit.nanometer
complex_settings.engine_settings.compute_platform = None

solvent_protocol = RelativeHybridTopologyProtocol(settings=solvent_settings)
complex_protocol = RelativeHybridTopologyProtocol(settings=complex_settings)
```

Verified defaults at OpenFE 1.13/main: `protocol_repeats=3`, equilibration 1 ns,
production 5 ns, `solvent_padding` 1.5 nm,
`engine_settings.compute_platform="cuda"`,
`alchemical_settings.softcore_LJ="gapsys"`.

- **The 10 ps / 50 ps values are a notebook cost hack.** Defaults are 1 ns
  equilibration and 5 ns production per repeat. If the user wants publishable
  numbers, agree on settings and convergence diagnostics with the user; defaults
  alone do not establish convergence. Never silently inherit demo lengths.
- **1.5 nm padding is meant for the solvent leg.** The notebook builds a second
  protocol with 1 nm padding for the complex leg, because 1.5 nm around a
  protein-ligand complex adds more water. These are tutorial choices, not a
  universal recommendation. Keep separate settings objects and record both paddings.
- `compute_platform` defaults to `"cuda"`; `None` lets OpenMMTools pick the best
  available platform, which is what you want without a CUDA GPU.
- OpenFE ≥1.13 also has a **private** `_adaptive_settings(stateA, stateB, mapping)`
  helper that widens the protocol for net-charge changes and narrows padding for
  protein systems. It is not public API and not used by the showcase — do not
  call it blindly.

### 6. Transformations and dry run

```python
import openfe

transformation_complex = openfe.Transformation(
    stateA=a_complex, stateB=b_complex, mapping=edge,
    protocol=complex_protocol, name=f"{a_complex.name}_{b_complex.name}_complex")
transformation_solvent = openfe.Transformation(
    stateA=a_solvent, stateB=b_solvent, mapping=edge,
    protocol=solvent_protocol, name=f"{a_solvent.name}_{b_solvent.name}_solvent")

complex_dag = transformation_complex.create()   # gufe ProtocolDAG
solvent_dag = transformation_solvent.create()

complex_unit = list(complex_dag.protocol_units)[0]
print(complex_unit.run(dry=True, verbose=True))  # {} on success, no sampling started
```

`Transformation.create()` builds the DAG of `ProtocolUnit`s; for this protocol the
units are the independent repeats. The dry run builds the hybrid OpenMM system and
stops before sampling — always do it before a long run.

### 7. Hand off to the CLI

```python
import pathlib
out = pathlib.Path("tyk2_json"); out.mkdir(exist_ok=True)
transformation_complex.to_json(out / f"{transformation_complex.name}.json")
transformation_solvent.to_json(out / f"{transformation_solvent.name}.json")
```

```bash
openfe quickrun tyk2_json/lig_ejm_31_lig_ejm_47_complex.json \
    -o results/lig_ejm_31_lig_ejm_47_complex.json \
    -d results/lig_ejm_31_lig_ejm_47_complex
```

This is the recommended production path: serialise, then run each leg under the
scheduler. See `openfe-rbfe-cli`.

### 8. Analysis after both legs have completed

Run the solvent transformation as well as the complex transformation, with matching
settings and independent repeats, before reporting an RBFE. The complex-leg energy
alone is not a relative binding free energy. Analyze existing results; do not start
`execute_DAG` under the guise of analysis (it launches simulations).

```bash
openfe gather results/ --report raw -o leg_results.tsv
openfe gather results/ --report ddg -o relative_results.tsv
```

Use `--report dg` only for a connected network with sufficient repeats. Its inferred
ligand values have an arbitrary additive offset, not a physical absolute affinity.
For Python analysis, call `complex_protocol.gather([completed_complex_dag_result])`
and `solvent_protocol.gather([completed_solvent_dag_result])` on results from existing,
approved runs. The notebook's `rbfe_protocol` variable is undefined; use the actual
protocol objects. Subtract solvent-leg estimates from complex-leg estimates with
the same A-to-B orientation, preserving units and propagating uncertainty.

## Validation and expected outputs

- Dry run returns `{}` with no traceback.
- `plot_atommapping_network` returns a `matplotlib.figure.Figure`; the graph has
  one node per ligand (10 for the tutorial set). GraphML round trip preserves
  nodes, edges, and mappings.
- Cheap pre-flight check: assert `len(dag.protocol_units) == protocol_repeats`
  for every transformation before spending GPU time.
- Stored tutorial `--report dg` values (not values from the short smoke-test\n  settings above, and not physical absolute affinities):

```text
ligand	DG(MLE) (kcal/mol)	uncertainty (kcal/mol)
lig_ejm_31	-0.09	0.05
lig_ejm_42	0.7	0.1
lig_ejm_46	-0.98	0.05
```

## Failure modes and safety

- **`StopIteration` from `suggest_mappings`** — no mapping for that pair with
  those mapper settings. Try the other mapper, or relax
  `map_exact_ring_matches_only` / `allow_bond_breaks` explicitly and say so.
- **`view_3d()` ImportError** — `py3Dmol` not installed; it is optional.
- **ImportError for `generate_lomap_network`** — the external `lomap2` package is
  missing.
- **Platform errors** — with `compute_platform="cuda"` on a CUDA-less machine
  OpenMM fails to create a context. Use `None`, `"OpenCL"`, or `"CPU"`.
- **Dry-run failure on lengths/conformers** — usually a ligand with no 3D
  conformer or missing partial charges. Check inputs before blaming the protocol.
- **`execute_DAG` is serial and in-process** — it blocks and holds memory for the
  whole DAG. Prototyping tool, not a scheduler; for real work serialise and
  dispatch `openfe quickrun`.
- **Cost safety** — hybrid-topology construction plus sampling is expensive.
  Confirm with the user before any non-dry-run simulation, before raising
  `protocol_repeats`, before raising production length, and before installing the
  OpenFEforge/condacolab environment in Colab.
- **Do not silently modernise scientific settings.** Where the notebook's values
  differ from current defaults (production length, Colab installer), report the
  difference and let the user choose.

## Sources

Showcase (rendered docs, checked 2026-10-26):
<https://docs.openfree.energy/en/latest/tutorials/showcase_notebook.html>

Upstream notebook, OpenFreeEnergy/ExampleNotebooks @
`7c1818c6d8954daf3df902fa6e5e224226dc0a5d` (main, 2026-10-26):
- <https://github.com/OpenFreeEnergy/ExampleNotebooks/blob/7c1818c6d8954daf3df902fa6e5e224226dc0a5d/showcase/openfe_showcase.ipynb>
- <https://github.com/OpenFreeEnergy/ExampleNotebooks/blob/7c1818c6d8954daf3df902fa6e5e224226dc0a5d/environment.yaml>
- Directory: <https://github.com/OpenFreeEnergy/ExampleNotebooks/tree/main/showcase>
  (inputs `showcase/inputs/tyk2_ligands.sdf` and `tyk2_protein.pdb` are
  byte-identical to the `rbfe_tutorial/` copies, md5 verified)

API verified against OpenFreeEnergy/openfe @
`16a5feaa0a7e19d2348d9d63ee15ea51fd471e15` (main, 2026-10-26); release `v1.13.0`:
`src/openfe/protocols/openmm_rfe/hybridtop_protocols.py`,
`src/openfe/protocols/openmm_rfe/equil_rfe_settings.py`,
`src/openfe/protocols/openmm_utils/omm_settings.py`,
`src/openfe/setup/ligand_network_planning.py`,
`src/openfe/utils/atommapping_network_plotting.py` —
<https://github.com/OpenFreeEnergy/openfe/tree/main/src/openfe>

`LigandAtomMapping.draw_to_file` / `view_3d`, `Protocol.create` / `Protocol.gather`,
`execute_DAG`, `ChemicalSystem`, `SolventComponent`, `ProteinComponent.from_pdb_file`,
`LigandNetwork.to_graphml/from_graphml` read from OpenFreeEnergy/gufe @
`996f92d81640cd9ca91beb8cbae7109f88fefb92` (main, 2026-10-26):
<https://github.com/OpenFreeEnergy/gufe/tree/main/src/gufe>

Install docs: <https://docs.openfree.energy/en/latest/installation.html>

Licensing: OpenFE, gufe, and the ExampleNotebooks are MIT; this skill is an
original summary, not a copy of upstream text.

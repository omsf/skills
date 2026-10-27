---
name: openadmet-data-catalogs
description: Discover, load, and correctly cite OpenADMET's curated ADMET datasets through Intake catalogs (ChEMBL, BindingDB, EveBio activities; Caco2/MDCK permeability, LogD, microsomal stability, PPB; tox21; RCSB/co-fold structures). Use when a user wants ADMET training data, asks how to read an OpenADMET catalog, hits a catalog-loading error, needs the raw-vs-aggregated schema, or needs per-source license and curation caveats before using the numbers.
---

# OpenADMET Data Catalogs

## Purpose and scope

`OpenADMET/data-catalogs` publishes curated ADMET datasets plus Intake catalog
YAML files that describe where each table lives. Catalogs cover:

- **Activity / ADMET endpoints** under `catalogs/activities/`: ChEMBL pChEMBL
  (Ki / IC50 / EC50 / permissive), ChEMBL Caco2 and MDCK permeability, ChEMBL
  LogD, ChEMBL microsomal stability, ChEMBL PPB, ChEMBL aminergic GPCR panel,
  BindingDB (Ki / IC50 / EC50 / permissive), EveBio (IC50 / EC50), and a
  PubChem tox21 AHR/PXR set.
- **Structures** under `catalogs/structures/`: RCSB mmCIF/PDB files and Boltz
  co-folding inputs. These are raw files — there are **no Intake YAML catalogs**
  for structures, so do not go looking for one.

Scope: obtaining and reading the data correctly. This skill does not cover
training models (see the OpenADMET demos/`openadmet-toolkit` for that).

Note up front: the repo README states it is "under very active development, we
make no guarantees about the stability or correctness of any catalogs contained
herein." Treat schemas as moving targets and verify at read time.

## Required inputs and clarifications

- **Which endpoint and which source?** e.g. "ChEMBL37 IC50 for PXR" vs
  "BindingDB permissive for LXR_alpha". Directory layout is
  `catalogs/activities/<SOURCE>_<ENDPOINT>/CATALOG_<VERSION>.yaml`.
- **raw or aggregated?** `*_raw.parquet` is one row per source measurement;
  `*_aggregated.parquet` is one row per canonicalised compound. Some families
  also ship `*_trimmed_aggregated`.
- **License acceptability.** Sources differ: ChEMBL is CC BY-SA 3.0, BindingDB
  CC BY 3.0, and **EveBio is CC BY-NC-SA 4.0 (non-commercial)**. Confirm with
  the user before using EveBio data in a commercial context.
- **Network/disk budget.** Reading a catalog pulls Parquet files from S3 on
  demand; the largest tables run to tens of thousands of rows. Cloning the repo
  copies ~2.4k files. Ask before bulk-downloading everything.

## Environment

Verified working with Python 3.13, `intake==2.0.9`, `pyarrow==26.0.0`,
`pandas==3.1.0`:

```bash
python -m venv .venv && source .venv/bin/activate
pip install intake pyarrow pandas aiohttp
```

Why each is needed (each verified empirically, in a fresh venv):

- `intake` — provides `intake.open_catalog` in 2.x; the catalogs are Intake
  **version 2** files (`version: 2`, `intake.readers.datatypes:Parquet`).
- `aiohttp` (+ `requests`) — fsspec's `HTTPFileSystem`. **Hard requirement only
  for remote catalogs.** Without it, `intake.open_catalog("https://...")` raises
  `ImportError: HTTPFileSystem requires "requests" and "aiohttp" to be installed`.
  Opening the same YAML from a local clone works without it.
- `pyarrow` / `pandas` — the reader is `intake.readers.readers:PandasParquet`
  with `output_instance: pandas:DataFrame`.

**`s3fs` is NOT required**, even though every data URL is
`s3://openadmet-data-public-dev/...`. pyarrow's built-in S3 filesystem handles
those URLs: a fresh environment with only `intake pyarrow pandas aiohttp`
(no `s3fs`, no `requests`) read both `CAR_aggregated` (226x15) and
`LXR_alpha_raw` (1720x30) correctly. Installing `s3fs` also works, but do not
add it as a documented prerequisite.

Note: the repo README says only `pip install intake`. That is not enough for the
URL-based workflow it demonstrates — a bare `intake` install fails with the
`HTTPFileSystem` ImportError above. Use the four-package install here.

No AWS credentials are needed. The objects are anonymously readable
(empirically checked for the example entries; `storage_options: null` does not
itself guarantee anonymous access); only bucket *listing* is
denied, which is irrelevant here.

## Procedure

### 1. Pick a catalog and build a raw URL

Browse the layout:

```bash
curl -s https://api.github.com/repos/OpenADMET/data-catalogs/contents/catalogs/activities \
  | python3 -c 'import json,sys; print("\n".join(e["name"] for e in json.load(sys.stdin)))'
```

Then list the YAML files and choose one:

```text
catalogs/activities/ChEMBL_pChEMBL_permissive/CATALOG_ChEMBL37_permissive.yaml
```

**Use `raw.githubusercontent.com`, not `github.com/blob/`.** The README's
example uses the `blob/main` HTML URL, which makes intake fetch an HTML page and
fail with `ScannerError: mapping values are not allowed here`. Verified working
form (pin a commit for reproducibility):

```text
https://raw.githubusercontent.com/OpenADMET/data-catalogs/<COMMIT_SHA>/catalogs/activities/ChEMBL_pChEMBL_permissive/CATALOG_ChEMBL37_permissive.yaml
```

### 2. Open the catalog and list entries

```python
import intake

CATALOG_URL = (
    "https://raw.githubusercontent.com/OpenADMET/data-catalogs/"
    "ac2234e842297895890d688cbb2284c3747d244d/"
    "catalogs/activities/ChEMBL_pChEMBL_permissive/CATALOG_ChEMBL37_permissive.yaml"
)

cat = intake.open_catalog(CATALOG_URL)
print(list(cat.entries))
```

Expected (verified): 24 entries, alternating `<TARGET>_aggregated` /
`<TARGET>_raw`:

```text
['AHR_aggregated', 'AHR_raw', 'CAR_aggregated', 'CAR_raw', 'CYP1A2_aggregated',
 'CYP1A2_raw', 'CYP2C9_aggregated', 'CYP2C9_raw', 'CYP2D6_aggregated',
 'CYP2D6_raw', 'CYP2J2_aggregated', 'CYP2J2_raw', 'CYP3A4_aggregated',
 'CYP3A4_raw', 'FXR_aggregated', 'FXR_raw', 'HERG_aggregated', 'HERG_raw',
 'LXR_alpha_aggregated', 'LXR_alpha_raw', 'PXR_aggregated', 'PXR_raw',
 'VDR_aggregated', 'VDR_raw']
```

`intake.open_catalog(...)` returns an `intake.readers.entry.Catalog`. Iterate
`cat.entries` (a mapping) — `list(cat)` returns `[]` and is a trap.

Local-file form works identically if you cloned the repo:

```python
cat = intake.open_catalog("data-catalogs/catalogs/activities/ChEMBL_pChEMBL_permissive/CATALOG_ChEMBL37_permissive.yaml")
```

### 3. Inspect an entry, then read it

```python
print(cat["PXR_aggregated"])     # reader class + kwargs + output instance + metadata
df = cat["PXR_aggregated"].read()
print(df.shape)                  # verified: (952, 9)
print(list(df.columns))
```

### 4. Verify the schema — do not trust the README

The READMEs document a single aggregated schema (`n_records` plus
`OPENADMET_TARGET`, `OPENADMET_TARGET_ID`, `OPENADMET_ENDPOINT`,
`OPENADMET_UNITS`, `OPENADMET_SOURCE`, `OPENADMET_SOURCE_VERSION`). **Two
different aggregated schemas coexist in the same catalog**, and this is easy to
get wrong. Verified:

| Entry | Rows | Cols | Schema |
| --- | --- | --- | --- |
| `PXR_aggregated` (ChEMBL37) | 952 | 9 | `OPENADMET_CANONICAL_SMILES`, `OPENADMET_INCHIKEY`, `assay_id_count`, `standard_value_{mean,median,std}`, `pchembl_value_{mean,median,std}` |
| `CYP1A2_aggregated` (ChEMBL37) | 11143 | 9 | same 9-column schema |
| `CAR_aggregated` (ChEMBL37) | 226 | 15 | `... n_records`, `OPENADMET_TARGET`, `OPENADMET_TARGET_ID`, `OPENADMET_ENDPOINT`, `OPENADMET_UNITS`, `OPENADMET_SOURCE`, `OPENADMET_SOURCE_VERSION` |
| `LXR_alpha_aggregated` (ChEMBL37) | 1207 | 15 | the documented 15-column schema |
| `LXR_alpha_raw` (ChEMBL37) | 1720 | 30 | documented raw schema |
| `VDR_raw` (ChEMBL37) | 847 | 30 | documented raw schema |

So: assert on the columns you need rather than selecting them blind.

```python
import pandas as pd

df = cat["LXR_alpha_aggregated"].read()
required = {"OPENADMET_CANONICAL_SMILES", "OPENADMET_INCHIKEY", "pchembl_value_mean"}
missing = required - set(df.columns)
if missing:
    raise RuntimeError(f"catalog schema drift: missing {missing}; got {list(df.columns)}")
```

Note the documented READMEs reflect the *newer* REST-API-built depositions
(CAR/FXR/VDR/LXR_alpha). Entries built by the older notebook path (AHR, PXR,
CYP1A2, CYP3A4, CYP2C9, CYP2D6, CYP2J2, HERG) still carry the 9-column schema.
Also note the raw schema is only described in the newer per-collection READMEs
(e.g. `ChEMBL_pChEMBL_permissive/README.md`), and `pchembl_value_std` /
`standard_value_std` are **null when `n_records == 1`** — that is expected, not
missing data.

### 5. Offline / bulk path

The repo commits the Parquet files alongside the YAMLs, so a clone gives you the
bytes without S3:

```bash
git clone --depth 1 https://github.com/OpenADMET/data-catalogs.git
```

```python
import pandas as pd

df = pd.read_parquet(
    "data-catalogs/catalogs/activities/ChEMBL_pChEMBL_permissive/"
    "ChEMBL37_permissive/ChEMBL_permissive_LXR_alpha_CHEMBL2808_aggregated.parquet"
)
print(df.shape)   # verified: (1207, 15), identical to the S3 read
```

Caveat: the committed catalog YAMLs contain **absolute `s3://` URLs**, so a
clone alone does *not* make `intake.open_catalog` offline — the README's claim
that you can "use the data directly by downloading it or cloning the repo"
applies to the Parquet files, not to the catalog YAMLs as shipped. Either read
the Parquet files directly (above) or rewrite the URLs to local paths.

### 6. Attribution

Quote the per-collection README for the exact attribution string and the
`OPENADMET_SOURCE_VERSION` value, for example:

```text
ChEMBL data is distributed under CC BY-SA 3.0 Unported.
Attribute as: ChEMBL (EMBL-EBI), release 37.
```

The repository itself is Apache-2.0; the *data* carries its own upstream license.
The archive is citable via DOI `10.5281/zenodo.21504232`.

## Validation and expected outputs

- `intake.open_catalog(...)` returns an `intake.readers.entry.Catalog`; a
  `ScannerError` means you fetched an HTML page (blob URL) or a truncated file.
- `len(list(cat.entries))` matches the count in the table below.
- `.read()` returns a `pandas.DataFrame` whose columns pass your explicit
  assertion from step 4.
- Cross-source sanity check: for a compound present in both ChEMBL and
  BindingDB, compare only matched targets, endpoints, units, assay contexts and
  censoring. In nM, an individual pActivity is `-log10(value * 1e-9)`, but the mean
  of logs is not the log of a mean: do not equate `pchembl_value_mean` with the
  transform of `standard_value_mean`. BindingDB's `standard_value*` are documented
  as nM, not specified by an API units field.

Catalog inventory at the verified revision (27 YAMLs, `entries` / `aggregated`
counts): BindingDB_EC50 8/4, BindingDB_IC50 8/4, BindingDB_Ki 2/1,
BindingDB_permissive 8/4, CHEMBL_CACO2 (35) 8/5, (37) 8/5, CHEMBL_MDCK (35) 6/3,
(37) 6/3, ChEMBL_GPCR_aminergic EC50 70/35, IC50 72/36, Ki 72/36, ChEMBL_LogD
(35) 3/2, (37) 3/2, ChEMBL_Microsomal (35) 6/3, (37) 6/3, ChEMBL_PPB (35) 6/4,
(37) 6/4, ChEMBL_pChEMBL_EC50 (35) 4/2, (37) 12/6, ChEMBL_pChEMBL_IC50 (35)
16/8, (37) 24/12, ChEMBL_pChEMBL_Ki (35) 16/8, (37) 18/9,
ChEMBL_pChEMBL_permissive (35) 16/8, (37) 24/12, EveBio_EC50 8/4,
EveBio_IC50 8/4.

## Failure modes and safety

- **`Mapping values are not allowed here` / `ScannerError`** — you used a
  `github.com/.../blob/...` URL. Switch to `raw.githubusercontent.com`.
- **`HTTPFileSystem requires "requests" and "aiohttp"`** — install `aiohttp`
  (only needed for remote catalogs).
- **`ImportError` / `No filesystem for protocol 's3'`** — unusual here; pyarrow
  reads `s3://` natively. If you hit it, install `s3fs` as a fallback rather
  than assuming it is required.
- **403 `AccessDenied` on the bucket root** — expected. Bucket listing is
  private; individual objects are readable. A 200 on a specific
  `https://openadmet-data-public-dev.s3.amazonaws.com/<key>` confirms access.
- **Column-mismatch between README and data** — see step 4. Assert explicitly.
- **ChEMBL version provenance is unreliable for the older notebook-built
  depositions.** The per-collection READMEs document that
  `curate_data_*.ipynb` calls `PermissiveChEMBLTargetCurator(version=...)`, but
  the real field is `chembl_version` (default `34`), and pydantic silently drops
  the unknown kwarg — so `ChEMBL35_*` / `ChEMBL37_*` directories built that way
  may actually contain ChEMBL34 data. The REST-API-built depositions
  (CAR/FXR/VDR/LXR_alpha in `ChEMBL37_*`, plus the BindingDB and EveBio sets) are
  not affected. Read `OPENADMET_SOURCE_VERSION` and the README before publishing.
- **Censored measurements are mixed into aggregates.** Aggregation ignores
  `standard_relation`, so `> 10000` and `< 1` records are averaged with exact
  ones. For compounds with many censored rows, go back to the `*_raw.parquet`
  `standard_relation` column.
- **EveBio** is CC BY-NC-SA 4.0 (non-commercial) and was retrieved from a
  **gated** Hugging Face dataset (`HF_TOKEN`); do not redistribute it casually,
  and do not use it in commercial work.
- **tox21** (`catalogs/activities/pubchem_tox21/`) has no Intake YAML — it ships
  CSV inputs, per-AID processed Parquet, `final/cleaned_*.parquet`, and
  `final/SCHEMA.md`. Read the Parquet with pandas.
- **`*_trimmed_aggregated`** entries (Caco2, LogD, PPB) are filtered variants of
  the full aggregate; check which one you actually meant.
- **Do not bulk-download the whole catalog set or clone full history** without
  asking; that is a real network/disk cost.

## Sources

Repository: <https://github.com/OpenADMET/data-catalogs>
Directory listing: <https://github.com/OpenADMET/data-catalogs/tree/main/catalogs/activities>

Verified at revision `ac2234e842297895890d688cbb2284c3747d244d`
(main, 2026-10-23, `Merge pull request #42 from OpenADMET/NHR-notebook-fix`):

- <https://github.com/OpenADMET/data-catalogs/blob/ac2234e842297895890d688cbb2284c3747d244d/README.md>
- <https://github.com/OpenADMET/data-catalogs/blob/ac2234e842297895890d688cbb2284c3747d244d/catalogs/activities/ChEMBL_pChEMBL_permissive/CATALOG_ChEMBL37_permissive.yaml>
- <https://github.com/OpenADMET/data-catalogs/blob/ac2234e842297895890d688cbb2284c3747d244d/catalogs/activities/ChEMBL_pChEMBL_permissive/README.md>
- <https://github.com/OpenADMET/data-catalogs/blob/ac2234e842297895890d688cbb2284c3747d244d/catalogs/activities/ChEMBL_pChEMBL_Ki/README.md>
- <https://github.com/OpenADMET/data-catalogs/blob/ac2234e842297895890d688cbb2284c3747d244d/catalogs/activities/EveBio_IC50/README.md>
- <https://github.com/OpenADMET/data-catalogs/blob/ac2234e842297895890d688cbb2284c3747d244d/catalogs/activities/BindingDB_IC50/README.md>
- <https://github.com/OpenADMET/data-catalogs/blob/ac2234e842297895890d688cbb2284c3747d244d/catalogs/activities/pubchem_tox21/final/SCHEMA.md>
- <https://github.com/OpenADMET/data-catalogs/blob/ac2234e842297895890d688cbb2284c3747d244d/catalogs/structures/README.md>
- <https://github.com/OpenADMET/data-catalogs/blob/ac2234e842297895890d688cbb2284c3747d244d/LICENSE> (Apache-2.0, repo)

Empirical checks performed while writing this skill (Python 3.13, intake 2.0.9,
pyarrow 26.0.0, pandas 3.1.0), each in a throwaway venv outside this repo:
all 27 catalog YAMLs open locally; remote open via `raw.githubusercontent.com`
works; `github.com/blob/` URL fails with `ScannerError`; `aiohttp` is a hard
requirement for remote catalogs and is absent-but-fine for local ones; `s3fs` is
not required (pyarrow reads the `s3://` Parquet directly); anonymous S3 object
GET returns 200 while bucket listing returns 403; S3 and local-clone row counts
and columns agree (1207x15); the 9-column vs 15-column aggregated split was
reproduced for PXR/CYP1A2 (9) and CAR/LXR_alpha (15).

Upstream data licenses are attributed in the per-collection READMEs above
(ChEMBL CC BY-SA 3.0, BindingDB CC BY 3.0, EveBio CC BY-NC-SA 4.0). This skill is
an original summary, not a copy of upstream text.
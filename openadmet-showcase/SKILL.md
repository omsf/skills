---
name: openadmet-showcase
description: End-to-end ADMET machine-learning workflow with the OpenADMET stack (openadmet-models / openadmet-toolkit) - curate ChEMBL activity data, train models with the Anvil YAML recipe system, compare models with post-hoc statistics, and predict on new compounds. Use for requests about the OpenADMET Showcase notebook, CYP3A4 inhibition modelling, Anvil recipes, `openadmet anvil|compare|predict`, or OpenADMET data catalogs.
---

# OpenADMET Showcase (Anvil training → comparison → inference)

## Purpose and scope

Reproduce and adapt the OpenADMET showcase pipeline: curate an ADMET endpoint
dataset (case study: CYP3A4 inhibition IC50 from ChEMBL), train candidate models
with the Anvil YAML recipe system, compare them with post-hoc statistics, and run
inference on unseen compounds. Scope covers the CLI workflows and the recipe format;
it does not cover force-field parameterization or free-energy calculations.

The upstream "showcase" notebook is a condensed version of the six numbered demos in
`openadmet-demos`; both are supported entry points.

## Required inputs and clarifications

Ask before proceeding if any of these are unknown:

1. **Endpoint / target**: ChEMBL target ID (CYP3A4 = `CHEMBL340`) and activity type
   (`IC50`, `Ki`, `EC50`, `AC50`, `Potency` are the accepted values).
2. **Data source**: pre-curated parquet from `OpenADMET/data-catalogs`, or a local
   curation run. Curation downloads the full ChEMBL SQLite database — confirm the
   disk/network budget and get explicit approval before starting it.
3. **Compute**: the LightGBM example runs on CPU; the shipped ChemProp/PyTorch
   Lightning recipes request GPUs (Slurm examples use `--gres=gpu:1 --time=48:00:00`).
   **Do not launch GPU training or the ChEMBL download without explicit approval.**
4. **Label name / task column**: `OPENADMET_LOGAC50` in the showcase, or
   `OPENADMET_LOGAC50_cyp3a4` for the multitask model.
5. Whether the user wants a single train/test split (fast, showcase default) or
   repeated k-fold cross validation (*required* for `openadmet compare`).

## Environment setup and version caveats

```bash
git clone https://github.com/OpenADMET/openadmet-demos.git
cd openadmet-demos
conda env create -f environment.yml      # env name is `openadmet-demos`
conda activate openadmet-demos
openadmet --help                         # verifies the CLI entry point
mkdir ../openadmet-showcase-run && cd ../openadmet-showcase-run
cp ../openadmet-demos/showcase/compounds_for_inference.csv .
```

Verified against upstream sources:

- `environment.yml` sets `name: openadmet-demos` and requires `python >= 3.12`. The
  README's `conda env create -n openadmet` followed by `conda activate openadmet-demos`
  is inconsistent; use `-f environment.yml`.
- Several critical packages come from git **`@main`** (openadmet-models,
  openadmet-toolkit, molfeat, mtenn, TabPFN, neural-pairwise-regression). Unpinned
  moving targets — record the resolved commit hashes for reproducibility.
- `zarr >=3.0.0` is required by the cached feature stores (`X_*_feat.zarr`); the data
  flow uses `intake` + `intake-parquet`.
- The Colab path uses the prebuilt installer
  `OpenADMETForge-0.1.0-alpha-Linux-x86_64.sh` from the `v0.1.0-alpha` release. The
  latest release tag `v0.2.0` has **no installer assets** at the time of writing.
- `openadmet` is a console script from `openadmet-models`; `docs/requirements.yaml` is
  only for building the Sphinx site.

## Steps

### 1. Obtain curated endpoint data

Pre-curated route (downloads only an endpoint Parquet file, not all of ChEMBL):

```python
import pandas as pd

url = (
    "https://raw.githubusercontent.com/OpenADMET/data-catalogs/"
    "ac2234e842297895890d688cbb2284c3747d244d"
    "/catalogs/activities/ChEMBL_pChEMBL_IC50/ChEMBL35_IC50"
    "/ChEMBL_IC50_CYP3A4_CHEMBL340_aggregated.parquet"
)
cyp3a4_data_agg = pd.read_parquet(url)
assert "OPENADMET_CANONICAL_SMILES" in cyp3a4_data_agg.columns
assert "pchembl_value_mean" in cyp3a4_data_agg.columns
```

Local curation route (downloads ChEMBL; **ask first**). Note the current
`PermissiveChEMBLTargetCurator` validator rejects the `require_units="uM"` spelling
shown in the showcase notebook markdown:

```python
from openadmet.toolkit.database.chembl import PermissiveChEMBLTargetCurator

curator = PermissiveChEMBLTargetCurator(
    chembl_version=35,          # or omit to auto-resolve the latest release
    chembl_target_id="CHEMBL340",
    standard_type="IC50",
    require_pchembl=True,
    require_units="nM",         # "nM" | "uM" | "mM"; must be "nM" with require_pchembl
)
cyp3a4_data_raw = curator.get_activity_data()
cyp3a4_data_agg = curator.aggregate_activity_data_by_compound()
```

### 2. Standardize the label columns and write the training CSV

```python
cyp3a4_data_agg["OPENADMET_ACTIVITY_TYPE"] = "IC50"
cyp3a4_data_agg["Target"] = "CYP3A4"
cyp3a4_data_agg["OPENADMET_LOGAC50"] = cyp3a4_data_agg["pchembl_value_mean"]
cyp3a4_data_agg.to_csv("cleaned_cyp3a4.csv", index=False)
```

### 3. Write an Anvil recipe

Minimum viable CPU recipe; use `<PLACEHOLDER>` values rather than inventing paths:

```yaml
data:
  resource: cleaned_cyp3a4.csv
  type: intake
  input_col: OPENADMET_CANONICAL_SMILES
  target_cols:
    - OPENADMET_LOGAC50
  dropna: true

metadata:
  authors: <YOUR NAME>
  email: <YOUR EMAIL>
  biotargets: [CYP3A4]
  build_number: 0
  description: basic regression using a LightGBM model
  driver: sklearn
  name: lgbm_pchembl
  tag: openadmet-chembl              # becomes the prediction column prefix
  tags: [openadmet, test, pchembl]
  version: v1

procedure:
  feat:
    type: FingerprintFeaturizer
    params:
      fp_type: "ecfp:4"
  model:
    type: LGBMRegressorModel
    params:
      alpha: 0.005
      learning_rate: 0.05
      n_estimators: 500
  split:
    type: ShuffleSplitter
    params:
      random_state: 42
      train_size: 0.8
      val_size: 0.0
      test_size: 0.2
  train:
    type: SKLearnBasicTrainer

report:
  eval:
    - type: RegressionMetrics
      params: {}
    - type: RegressionPlots
      params:
        axes_labels: [True pAC50, Predicted pAC50]
        max_val: 10
        min_val: 3
        pXC50: true
        title: True vs Predicted pAC50 on test set
```

Swap in `SKLearnRepeatedKFoldCrossValidation` (with `n_splits`/`n_repeats`) for the
report stage when the results feed `openadmet compare`. The demo's richer recipe uses
`FeatureConcatenator(DescriptorFeaturizer desc2d + FingerprintFeaturizer ecfp:4)` —
a modeling choice, not a requirement, and it changes runtime.

### 4. Train, compare, predict

```bash
export OADMET_NO_RICH_LOGGING=1         # friendlier output inside notebooks/CI
openadmet anvil --recipe-path anvil_lgbm.yaml --output-dir lgbm
```

Save the recipe above as `anvil_lgbm.yaml` first, replacing the name/email placeholders.
For the optional comparison below, train the CV recipe into `lgbm_cv` and train the
multitask model with matching folds. The single-split recipe above does not produce
these CV artifacts. Do not silently use the repository's precomputed models as if
they were the model you just trained.

```bash
export OADMET_NO_RICH_LOGGING=1
mkdir -p model_comparisons
openadmet compare \
    --model-stats-fns ./lgbm_cv/cross_validation_metrics.json \
    --labels lgbm \
    --task-names OPENADMET_LOGAC50 \
    --model-stats-fns ./multitask/cross_validation_metrics.json \
    --labels multitask \
    --task-names OPENADMET_LOGAC50_cyp3a4 \
    --mt-id cyp3a4 \
    --output-dir model_comparisons/ \
    --report True
```

The notebook markdown documents an older flag set (`--model-stats`, `--model-tag`,
`--task-name`); the shipped CLI uses `--model-stats-fns`, `--labels`, `--task-names`.
Prefer the full-directory form (`--model-dirs` + `--label-types`, with label types
`biotarget|model|feat|tasks`) when the whole Anvil output tree is available.

```bash
export OADMET_NO_RICH_LOGGING=1
openadmet predict \
    --input-path compounds_for_inference.csv \
    --input-col OPENADMET_CANONICAL_SMILES \
    --model-dir lgbm \
    --output-csv predictions.csv \
    --accelerator cpu
```

This predicts with the locally trained `lgbm` model. `--model-dir` is repeatable:
pass it once per trusted, locally trained ensemble member. Model artifacts may
contain pickle objects and execute code when loaded; never load untrusted models.

## Validation and expected outputs

Anvil output directory (`--output-dir`) contains `anvil_recipe.yaml`,
`recipe_components/{data,metadata,procedure,eval}.yaml`,
`data/{X_train,X_val,X_test,y_train,y_val,y_test}.csv` plus `*_feat.zarr`,
`model.json`, `model.pkl` (or `model.pth` for ChemProp), `regression_metrics.json`,
and — when the report section includes CV — `cross_validation_metrics.json` and
`cross_validation_regplot.png`.

`openadmet compare` writes `anova.pdf`, `Levene.json`, `Tukey_HSD.json`,
`mcs_plots.pdf`, `mean_diffs.pdf`, `normality_plots.pdf`, `paired_plots.pdf`,
`posthoc.pdf` into `--output-dir`, plus a PDF report when `--report True`.

`openadmet predict` appends `OADMET_PRED_<metadata.tag>_<task_column>` and
`OADMET_STD_<metadata.tag>_<task_column>` (std is NaN for a single model, populated for
ensembles), so the showcase output column is
`OADMET_PRED_openadmet-chembl_OPENADMET_LOGAC50`. Sanity checks: row count matches the
input, no NaNs in the PRED column, and the predicted `pIC50` distribution overlaps the
training distribution.

## Failure modes and safety

- **Comparing incomparable models**: all models must share identical CV structure
  (e.g. 5 splits × 5 repeats); mismatched folds make the statistics meaningless.
- **Multitask comparison**: `--mt-id` is required when compared models differ in task
  composition; choose an identifier absent from every other model's target columns.
- **`require_units` validator**: with `require_pchembl=True` only `"nM"` is accepted;
  anything else raises `ValueError` at construction.
- **Do not silently change `random_state`, split sizes, or CV folds** when reproducing
  published numbers — they change the reported metrics.
- **GPU cost**: ChemProp/multitask recipes set `accelerator: gpu`, `max_epochs: 50` with
  early stopping. Confirm with the user before running; use `--accelerator cpu`
  inference for small screens.
- **CPU inference hang**: the finetuning demo reports the inference CLI hanging on CPU
  for deep-learning models *after* writing the CSV. Do not infer success from file
  existence alone: verify row count, finite predictions, and error logs before asking
  the user whether to interrupt a hung process.
- **Private S3 intake catalogs**: fall back to the GitHub-hosted parquet files.
- **ChEMBL version drift**: the showcase pins ChEMBL 35; newer releases (e.g. 37) exist
  in `data-catalogs`. Changing the release changes the dataset — record it.

## Sources

Read at commit `8b2c5f34ac892a3aa256d8564bd33c707fe9d11b`
(`OpenADMET/openadmet-demos`, 2026-07-23), directory
<https://github.com/OpenADMET/openadmet-demos>:

- Showcase notebook: `showcase/OpenADMET_Models_Showcase.ipynb`
  <https://github.com/OpenADMET/openadmet-demos/blob/8b2c5f34ac892a3aa256d8564bd33c707fe9d11b/showcase/OpenADMET_Models_Showcase.ipynb>
- Recipe: `showcase/anvil_lgbm.yaml`; demos under `demos/01..06`
- Env/CI: `environment.yml`, `docs/requirements.yaml`,
  `.github/workflows/ci.yaml` (runs notebooks via `pytest --nbmake`),
  `.github/workflows/make_installer.yaml`

CLI/API verified against `OpenADMET/openadmet-models` @
`30d491bf55eabc2e8a1aefa6d49af411ef179280` (2026-10-23):
`openadmet/models/cli/{cli,anvil,compare,predict}.py`,
`openadmet/models/comparison/posthoc.py`,
`openadmet/models/inference/inference.py`; and `OpenADMET/openadmet-toolkit` @
`a7c110e6dde468afd9a3a3d84bab6d3be2e45467` (2026-10-14):
`openadmet/toolkit/database/chembl.py`.

Data: `OpenADMET/data-catalogs` @ `ac2234e842297895890d688cbb2284c3747d244d`
(2026-10-23), `catalogs/activities/ChEMBL_pChEMBL_IC50/`.

Licenses: Apache-2.0 (`openadmet-demos`, `data-catalogs`, verified in `LICENSE`). Cite
upstream projects and the Zenodo DOIs advertised in their READMEs when reusing data.

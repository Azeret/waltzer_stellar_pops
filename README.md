# waltzer_stellar_pops

Helper scripts + example outputs to build **WALTzER ETC** input target lists for nearby stellar populations (open clusters) and summarize feasibility.

This repo **does not** bundle the ETC itself; it assumes you also have the WALTzER ETC repo available.

## Expected folder layout

Recommended: clone this repo next to the ETC:

```text
workdir/
  waltzer_stellar_pops/   (this repo)
  waltzer_etc/            (https://github.com/pcubillos/waltzer_etc)
  visibility/             (optional, https://github.com/AlexisBrandeker/visibility)
```

## Quick start

### 1) Python environment (for target-list building)

```bash
python -m venv .venv
.venv/bin/python -m pip install -U pip
.venv/bin/python -m pip install -r requirements.txt
```

### 2) Build a cluster target list (Gaia DR3)

Example: M45 (Pleiades) + NGC 2632 (Praesepe) + ONC:

```bash
.venv/bin/python scripts/build_cluster_targets.py \
  --cluster "M45" --cluster "NGC 2632" --cluster "ONC" \
  --radius-deg 1.5 \
  --g-max 12.5 \
  --obs-hours 2.5 \
  --out targets_open_clusters.csv
```

Output is a `waltz`-compatible CSV with required columns:
`pl_name, pl_trandur, st_teff, st_rad, st_mass, ra, dec, sy_vmag`
plus extra diagnostic columns (Gaia IDs, parallax/PM, neighbor metrics, and a Sun-angle-only visibility window).

## Running the WALTzER ETC

### 1) Install the ETC

From the sibling `waltzer_etc/` repo:

```bash
cd ../waltzer_etc
python -m venv .venv
.venv/bin/python -m pip install -U pip
.venv/bin/python -m pip install -e .
```

### 2) Download the stellar SED libraries

The ETC supports `--sed llmodels`, `--sed phoenix`, and `--sed bt_settl`, but the grids are not in the repo by default.

This repo includes an installer script that downloads and places templates in the folders expected by WALTzER:

```bash
cd ../waltzer_stellar_pops
scripts/download_sed_templates.sh --waltzer-etc-dir ../waltzer_etc --llmodels
```

Notes:
- PHOENIX is a **large** download (~1.8 GB).
- LLMODELS is downloaded from a Google Drive link used by the WALTzER ETC docs.

### 3) Run an ETC feasibility pass

Example (2.5h in-transit, 1 visit, 35 cm telescope, PHOENIX SEDs):

```bash
cd ../waltzer_stellar_pops
../waltzer_etc/.venv/bin/waltz targets_open_clusters.csv waltzer_snr.csv \
  --sed phoenix --tdur 2.5 --nobs 1 --diam 35 --eff 0.6
```

## Examples and docs

- Example target list: `examples/targets_open_clusters_2026-02-06.csv`
- Example ETC results: `examples/waltzer_snr_phoenix_tdur2p5_nobs1_2026-02-06.csv`
- Example merged feasibility table: `examples/feasibility_phoenix_tdur2p5_nobs1_merged_2026-02-06.csv`
- Feasibility summary: `docs/feasibility_summary_phoenix_tdur2p5_nobs1_2026-02-06.md`
- Work summary: `docs/2026-02-06_stellar_science_work_summary.md`

## Proposal figures (stellar science)

The current proposal-level stellar-science feasibility figures (A4-ready) and the minimal scripts/data needed to reproduce them live in:
- `stellar_science/README.md`

## Caveats

- ONC extinction/reddening is **not** modeled in the quick ETC pass; NUV feasibility there is likely optimistic.
- Neighbor/crowding metrics (`nn_*`) are computed against the Gaia cone-search magnitude limit you choose (`--g-max`). If you want fainter contaminants, re-run with a deeper `--g-max` and then down-select by brightness.

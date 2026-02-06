# Stellar science case: clusters + feasibility (2026-02-06)

This document summarizes the “nearby open cluster” feasibility workflow implemented in `waltzer_stellar_pops` and the outputs produced on **February 6, 2026**.

## Repositories used

- WALTzER ETC (CLI `waltz`): `pcubillos/waltzer_etc`
- Visibility checker (Sun-angle-only): `AlexisBrandeker/visibility` (optional; the target builder also computes the same style of visibility window internally)
- This repo: `Azeret/waltzer_stellar_pops`

## What was done

### 1) Built a cluster target list (Gaia DR3)

Script: `scripts/build_cluster_targets.py`

Clusters:
- M45 (Pleiades)
- NGC 2632 (Praesepe / Beehive)
- ONC (Orion Nebula Cluster)

Query + selection (first-pass, “simple”):
- Gaia DR3 cone search around the SIMBAD cluster center.
- `phot_g_mean_mag <= 12.5` (bright sample for quick feasibility).
- Require non-null Gaia values for `teff_gspphot`, `parallax`, `pmra`, `pmdec`.
- RUWE filter: `ruwe <= 1.4`.
- Membership cut: robust median/MAD window in parallax and proper motion (tolerances include floors so tiny MADs don’t over-tighten).
- Add nearest-neighbor metrics (`nn_sep_arcsec`, `nn_dG`) against *all* Gaia sources returned by the cone query (i.e., not just “members”).
- Compute an approximate Johnson V (`sy_vmag`) from Gaia `G` and `BP-RP` (polynomial transform; not extinction-corrected).
- Compute a Sun-angle-only visibility window string for year 2026 using a minimum solar separation of 110° (WALTzER-like).

Command used:

```bash
python scripts/build_cluster_targets.py \
  --cluster "M45" --cluster "NGC 2632" --cluster "ONC" \
  --radius-deg 1.5 \
  --g-max 12.5 \
  --obs-hours 2.5 \
  --out targets_open_clusters.csv
```

Output:
- `examples/targets_open_clusters_2026-02-06.csv` (983 targets total)

### 2) Downloaded/installed SED templates for the ETC

Script: `scripts/download_sed_templates.sh`

Installed:
- PHOENIX (`phoenixm00_*.fits`) into `../waltzer_etc/waltzer_etc/data/phoenix/`
- BT-Settl (`phoenixm0.0_*_5.0_2011.fits`) into `../waltzer_etc/waltzer_etc/data/bt_settl/`
- LLMODELS (`t*g4.4.flx`) into `../waltzer_etc/waltzer_etc/data/models/`

### 3) Ran the WALTzER ETC (quick feasibility pass)

ETC settings:
- `--sed phoenix`
- `--tdur 2.5` hours (fixed for all targets)
- `--nobs 1`
- `--diam 35` cm
- `--eff 0.6`

Outputs saved:
- `examples/waltzer_snr_phoenix_tdur2p5_nobs1_2026-02-06.csv`
- `docs/feasibility_summary_phoenix_tdur2p5_nobs1_2026-02-06.md`
- `examples/feasibility_phoenix_tdur2p5_nobs1_merged_2026-02-06.csv`

## Key numbers from the quick ETC run

All numbers below correspond to the PHOENIX run above (2.5h in-transit, 1 visit).

Targets per cluster (in the example list):
- M45: 401 targets; visibility window (2026, Sun-angle only): `0914-0130`
- NGC 2632: 380 targets; visibility window (2026, Sun-angle only): `1120-0407`
- ONC: 202 targets; visibility window (2026, Sun-angle only): `1010-0218`

“How many are good enough?” depends on the science thresholds, but as a quick sanity check:
- VIS precision: number of targets with `VIS_transit_uncert < 10,000 ppm`
  - M45: 229 / 401
  - NGC 2632: 212 / 380
  - ONC: 124 / 202
- NUV precision: number of targets with `NUV_transit_uncert < 3,000 ppm`
  - M45: 18 / 401
  - NGC 2632: 3 / 380
  - ONC: 12 / 202

## High-level results (interpretation)

- VIS: many targets in all three clusters reach high spectroscopic SNR for a 2.5h in-transit integration (per WALTzER sampling element).
- NIR: for this bright sample, SNR is very high across the board (photometric band-integrated in the ETC).
- NUV: feasibility is limited to the brightest/hottest stars; ONC has some very bright OB stars that look strong in NUV **in this simplified run**.

## Critical caveats (for sharing)

- ONC extinction/reddening is not modeled here. Because extinction rises strongly into the NUV, ONC NUV feasibility is likely **optimistic** unless you incorporate `A_V/E(B−V)` (either by dereddening inputs or reddening the SEDs used by the ETC).
- Crowding/blending risk is not modeled physically; `nn_sep_arcsec`/`nn_dG` is only a quick screen and depends on the Gaia magnitude limit used for the cone search (`--g-max`).

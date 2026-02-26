# Reproduce the proposal figures (stellar science)

This folder is intentionally kept minimal: it contains **only** the scripts + inputs needed to reproduce the figures currently used in the stellar-science theme.

## Figures in use

1) **ETC example spectra by category (A4)**
- `stellar_science/fig_spectra_by_category_variable_t_A4.png`
- `stellar_science/fig_spectra_by_category_variable_t_A4.pdf`

2) **All-sky Gaia DR3 density + reachability tiers (A4)**
- `stellar_science/fig_gaia_allsky_vdist_density_reach_snrproxy_A4.png`
- `stellar_science/fig_gaia_allsky_vdist_density_reach_snrproxy_A4.pdf`
3) **All-sky Gaia DR3 young-star density + reachability tiers (A4)**
- `stellar_science/fig_gaia_young_vdist_density_reach_snrproxy_A4.png`
- `stellar_science/fig_gaia_young_vdist_density_reach_snrproxy_A4.pdf`

## Inputs kept for reproducibility

**ETC-derived table + spectra**
- `stellar_science/data/etc/feasibility_phoenix_tdur2p5_nobs1_merged.csv`
  - contains Gaia-derived target metadata + ETC SNR summaries (baseline 2.5 h)
- `stellar_science/data/etc/waltzer_snr_phoenix_tdur2p5_nobs1.pickle`
  - per-wavelength spectra + variances used to render the NUV/VIS spectra panels

**Gaia histogram cache (optional but recommended)**
- `stellar_science/data/gaia_hist_snrproxy/`
  - binned-count query results used by the all-sky density plot
  - if missing, the script will re-query Gaia DR3 via `astroquery`

**Young-star membership input (Figure 3)**
- `stellar_science/data/young/banyan_members.tsv`
  - downloaded from VizieR (BANYAN Σ members; `J/ApJ/856/23/members`)
- `stellar_science/data/young/banyan_members_gaia_dr3_cache.csv`
  - optional cache produced by the script (crossmatch + Gaia photometry/Teff)

## Figure 1: ETC spectra by stellar-science category

Script:
- `stellar_science/make_feasibility_fig_spectra_categories.py`

Run:
```bash
python stellar_science/make_feasibility_fig_spectra_categories.py \
  --merged-csv stellar_science/data/etc/feasibility_phoenix_tdur2p5_nobs1_merged.csv \
  --pickle stellar_science/data/etc/waltzer_snr_phoenix_tdur2p5_nobs1.pickle \
  --out stellar_science/fig_spectra_by_category_variable_t_A4.png
```

Realistic wavelength spacing + full-SED overlay:
```bash
python stellar_science/make_feasibility_fig_spectra_categories.py \
  --merged-csv stellar_science/data/etc/feasibility_phoenix_tdur2p5_nobs1_merged.csv \
  --pickle stellar_science/data/etc/waltzer_snr_phoenix_tdur2p5_nobs1.pickle \
  --x-layout realistic \
  --out stellar_science/fig_spectra_by_category_variable_t_realx_A4.png
```

Assumptions encoded in the figure (proposal defaults):
- Baseline ETC exposure for SNR columns: **2.5 h**
- Variable-time scaling: **S/N ∝ √t**; per-target `t_need = max(NUV,VIS,NIR)` to meet the per-case SNR thresholds
- Max per-target time for “limit example” selection: **5 h**
- PHOENIX photospheric templates; normalization is by **V magnitude**
- NUV+VIS are **spectroscopy** (R≈3000); NIR is **band-integrated photometry** (single broad band)
- ONC example uses **no extinction** (optimistic for NUV)

Per-category SNR requirements used (median):
- Massive/hot: **NUV≥50, VIS≥100, NIR≥30**
- Quiet FGK: **NUV≥20, VIS≥200, NIR≥30**
- Active K/M: **NUV≥20, VIS≥100, NIR≥30**
- PMS K/M (ONC): **NUV≥5, VIS≥300, NIR≥30** (pragmatic; intended for rebinned/indices; no A_V in this demo)

## Figure 2: Gaia DR3 all-sky density + reachability tiers (≤1 h and ≤5 h)

Script:
- `stellar_science/gaia_allsky_vdist_density_reachability_snrproxy.py`

Run (uses cached histograms if present, otherwise queries Gaia DR3):
```bash
python stellar_science/gaia_allsky_vdist_density_reachability_snrproxy.py \
  --etc-merged stellar_science/data/etc/feasibility_phoenix_tdur2p5_nobs1_merged.csv \
  --cache-dir stellar_science/data/gaia_hist_snrproxy \
  --out-prefix stellar_science/fig_gaia_allsky_vdist_density_reach_snrproxy_A4
```

Optional: write a small markdown summary (kept out of the proposal folder by default):
```bash
python stellar_science/gaia_allsky_vdist_density_reachability_snrproxy.py \
  --etc-merged stellar_science/data/etc/feasibility_phoenix_tdur2p5_nobs1_merged.csv \
  --cache-dir stellar_science/data/gaia_hist_snrproxy \
  --out-prefix stellar_science/fig_gaia_allsky_vdist_density_reach_snrproxy_A4 \
  --out-md stellar_science/_archive/gaia_allsky_vdist_density_reachability_snrproxy.md
```

What “reachable in all bands” means in this figure:
- We fit an **empirical SNR proxy model** `SNR_band(V, Teff)` to the local ETC merged table (baseline 2.5 h).
- For each Gaia histogram bin `(distance, V_est, Teff)` we predict baseline SNRs in **NUV, VIS, NIR** and compute
  `t_need = max( t_need_NUV, t_need_VIS, t_need_NIR )` using **S/N ∝ √t**.

## Figure 3: Gaia DR3 young-star density + reachability tiers (≤1 h and ≤5 h)

Script:
- `stellar_science/gaia_allsky_young_vdist_density_reachability_snrproxy.py`

Run:
```bash
python stellar_science/gaia_allsky_young_vdist_density_reachability_snrproxy.py \
  --etc-merged stellar_science/data/etc/feasibility_phoenix_tdur2p5_nobs1_merged.csv \
  --out-prefix stellar_science/fig_gaia_young_vdist_density_reach_snrproxy_A4
```

Notes:
- “Young” is selected from the literature-based BANYAN Σ membership compilation (Gagné+ 2018; VizieR `J/ApJ/856/23/members`).
- The two panels show members in associations with approximate ages **≤100 Myr** and **≤10 Myr** (age mapping is encoded in the script).
- The script crossmatches BANYAN members to Gaia DR3 via the 2MASS best-neighbour table (`gaiadr3.tmass_psc_xsc_best_neighbour`).
- A small per-member Gaia cache is written to `stellar_science/data/young/banyan_members_gaia_dr3_cache.csv` to speed re-runs (safe to delete).
- Reachability tiers are computed the same way as Figure 2, using the ETC-trained SNR proxy and assuming **S/N ∝ √t**:
  - **≤1 h** (solid contour)
  - **≤5 h** (dashed contour)

### Gaia DR3 query (ADQL) used to build the all-sky histograms

The script queries `gaiadr3.gaia_source` via **Gaia TAP** (through `astroquery.gaia`) and returns only **binned counts**
in `(distance, V_est, Teff)`; it does **not** download per-source tables for these plots.

**V estimate used (no extinction correction):**
- `V_est = G + (0.01760 + 0.006860*(BP−RP) + 0.1732*(BP−RP)^2 − 0.045858*(BP−RP)^3)`
  - implemented in ADQL as `phot_g_mean_mag + (...)` using `bp_rp`.

**Histogram query template (defaults: Δd=10 pc, ΔV=0.25 mag, ΔTeff=200 K):**
```sql
SELECT dist_bin_pc, v_bin, teff_bin_k, COUNT(*) AS n
FROM (
  SELECT
    FLOOR((1000.0/parallax)/10.0)*10.0 AS dist_bin_pc,
    FLOOR((V_est)/0.25)*0.25 AS v_bin,
    FLOOR((teff_gspphot)/200.0)*200.0 AS teff_bin_k
  FROM gaiadr3.gaia_source
  WHERE <CASE_WHERE>
) AS sub
GROUP BY dist_bin_pc, v_bin, teff_bin_k
```

**Base quality cuts (defaults; see script flags to change):**
- `phot_g_mean_mag`, `bp_rp`, `teff_gspphot`, `logg_gspphot`, `parallax`, `ruwe`, `parallax_over_error` not null
- `parallax > 0`
- `phot_g_mean_mag ≤ 20` (guard)
- `ruwe ≤ 1.4`
- `parallax_over_error ≥ 10`
- `V_est ∈ [4,16]` (plot window)

**Per-science-case proxy selections (`<CASE_WHERE>`; defaults):**
- Quiet FGK dwarfs: `parallax ≥ 2.0 mas` (≤500 pc), `5200 ≤ teff_gspphot ≤ 6500`, `logg_gspphot ≥ 4.0`
- Active K/M dwarfs (activity/flares proxy): `parallax ≥ 2.0 mas` (≤500 pc), `teff_gspphot < 5200`, `logg_gspphot ≥ 4.0`
- PMS K/M candidates (proxy; not youth-selected): `parallax ≥ 2.0 mas` (≤500 pc), `teff_gspphot < 5200`, `3.2 ≤ logg_gspphot ≤ 4.0`
- Massive/hot: `parallax ≥ 0.5 mas` (≤2000 pc), `teff_gspphot ≥ 15000`

Important caveats:
- `V_est` is computed from Gaia **G and BP−RP** (polynomial), **no extinction correction**.
- The Gaia category selections use `teff_gspphot` and `logg_gspphot` and are **proxies**, not a vetted astrophysical classification.
- The SNR proxy is trained on the **current ETC sample in this repo**; treat it as proposal-level guidance.

“Massive/hot” definition used here:
- Gaia **`teff_gspphot ≥ 15000 K`** (proxy for O/B-type hot stars; not a strict spectral-type catalogue).

## Python dependencies (for reproducibility)

Both scripts use standard scientific Python:
- `numpy`, `pandas`, `matplotlib`

Figure 2 additionally needs:
- `astroquery` (only if you need to refresh the Gaia histogram cache)

## Exoplanet targets: population + visibility (from `waltzer_planets.xlsx`)

This produces two figures for the targets in `stellar_science/waltzer_planets.xlsx`:
- `stellar_science/fig_planets_population_snr_A4.png` (+ `.pdf`)
- `stellar_science/fig_planets_visibility_2026_A4.png` (+ `.pdf`)
- `stellar_science/fig_planets_spectra_by_category_variable_t_A4.png` (+ `.pdf`)
- `stellar_science/fig_planets_spectra_by_category_variable_t_realx_A4.png` (+ `.pdf`)

It also writes intermediate products to:
- `stellar_science/data/planets/`

Run (recommended via a local venv so `openpyxl` is available for `.xlsx`):
```bash
python -m venv --system-site-packages stellar_science/.venv
stellar_science/.venv/bin/python -m pip install openpyxl synphot pyratbay
stellar_science/.venv/bin/python stellar_science/make_planet_targets_population_visibility.py
stellar_science/.venv/bin/python stellar_science/make_planet_feasibility_fig_spectra_categories.py
```

Planet-host spectra with realistic wavelength spacing + full-SED overlay:
```bash
stellar_science/.venv/bin/python stellar_science/make_planet_feasibility_fig_spectra_categories.py \
  --x-layout realistic \
  --out stellar_science/fig_planets_spectra_by_category_variable_t_realx_A4.png
```

Notes:
- The script runs the ETC stage-1 using the local sibling repo `waltzer_etc/` (it is imported from source; not pip-installed).
- Visibility is Sun-angle only (default: year 2026, min Sun separation 110°).
- Population “reachability” is computed by tuning the number of stacked transits per target assuming `S/N ∝ √t`, using the same per-band SNR thresholds as the stellar-science plots (category-dependent).

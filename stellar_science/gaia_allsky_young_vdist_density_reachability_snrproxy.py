#!/usr/bin/env python3
"""
Gaia DR3 young stars: distance–V density with *SNR-based* reachability tiers (proxy).

Definition of "young" here
--------------------------
Gaia DR3 does not provide a robust all-sky ≤100 Myr selection via FLAME ages: the DR3
`age_flame` parameter effectively bottoms out around ~0.2 Gyr (200 Myr) for this use.

For a proposal-level, *astrophysically clean* youth selection we therefore use the
literature-based membership compilation from BANYAN Σ:
  - Gagné+ 2018, VizieR: J/ApJ/856/23/members

We assign an approximate age (Myr) per association and build two subsets:
  - ≤100 Myr
  - ≤10 Myr

We crossmatch members to Gaia DR3 via the 2MASS best-neighbour table, then compute:
  - V_est from Gaia (G, BP-RP; polynomial; no extinction correction)
  - distance from parallax
  - per-target time need from an ETC-trained SNR proxy model, assuming S/N ∝ √t

Important caveats
-----------------
This "young-star population" is *not complete*: it reflects known nearby associations
in the BANYAN Σ compilation (little/no embedded population; extinction not modelled).
Treat it as a clean lower bound / sanity check, not an all-sky census of ≤100 Myr stars.
"""

from __future__ import annotations

import argparse
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd


def _set_mpl_cache() -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(Path("stellar_science/.mplconfig").resolve()))


def _v_from_g_bp_rp(g: np.ndarray, bp_rp: np.ndarray) -> np.ndarray:
    g = np.asarray(g, dtype=float)
    c = np.asarray(bp_rp, dtype=float)
    return g + (0.01760 + 0.006860 * c + 0.1732 * c**2 - 0.045858 * c**3)


# Association ages (Myr): compact mapping for youth binning.
# Values are representative literature ages used widely in the nearby-YMG context.
ASSOC_AGE_MYR: dict[str, float] = {
    # Star-forming regions / very young (≲10 Myr)
    "ROPH": 1.0,  # ρ Oph
    "TAU": 1.5,  # Taurus
    "EPSC": 4.0,  # ε Cha
    "CRA": 4.5,  # CrA
    "USCO": 10.0,  # Upper Sco
    "TWA": 10.0,  # TW Hya
    "UCRA": 10.0,  # Upper CrA
    "ETAC": 11.0,  # η Cha
    # Sco-Cen older subgroups / young MGs (≲100 Myr)
    "LCC": 15.0,
    "UCL": 16.0,
    "THOR": 22.0,
    "BPMG": 24.0,
    "OCT": 35.0,
    "COL": 42.0,
    "CAR": 45.0,
    "THA": 45.0,
    "IC2391": 50.0,
    "IC2602": 46.0,
    # Groups present in the BANYAN table but older than our ≤100 Myr cut
    "ABDMG": 149.0,
    "PLE": 125.0,
    "HYA": 650.0,
    "UMA": 414.0,
    "CBER": 600.0,
    "CARN": 200.0,
    # Unknown/uncertain in this quick mapping (excluded by default)
    # "118TAU": ...
    # "PL8": ...
    # "XFOR": ...
}


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--etc-merged",
        default="stellar_science/data/etc/feasibility_phoenix_tdur2p5_nobs1_merged.csv",
        help="ETC merged table used to fit the SNR proxy model.",
    )
    p.add_argument(
        "--members-tsv",
        default="stellar_science/data/young/banyan_members.tsv",
        help="BANYAN Σ members TSV downloaded from VizieR (J/ApJ/856/23/members).",
    )
    p.add_argument(
        "--cache-members-gaia",
        default="stellar_science/data/young/banyan_members_gaia_dr3_cache.csv",
        help="Per-member Gaia DR3 cache (speeds re-runs; safe to delete).",
    )
    p.add_argument("--out-prefix", default="stellar_science/fig_gaia_young_vdist_density_reach_snrproxy_A4")
    p.add_argument("--dmax-pc", type=float, default=500.0)
    p.add_argument("--vmin-plot", type=float, default=4.0)
    p.add_argument("--vmax-plot", type=float, default=16.0)
    p.add_argument("--bin-dist-pc", type=float, default=10.0)
    p.add_argument("--bin-vmag", type=float, default=0.25)

    # Gaia quality cuts
    p.add_argument("--gmax-guard", type=float, default=20.0)
    p.add_argument("--ruwe-max", type=float, default=1.4)
    p.add_argument("--plxsnr-min", type=float, default=10.0)

    # SNR requirements (defaults to the "PMS" science-case thresholds used elsewhere in this repo)
    p.add_argument("--req-nuv", type=float, default=5.0)
    p.add_argument("--req-vis", type=float, default=300.0)
    p.add_argument("--req-nir", type=float, default=30.0)

    p.add_argument("--t1-hours", type=float, default=1.0)
    p.add_argument("--t5-hours", type=float, default=5.0)
    p.add_argument("--etc-base-hours", type=float, default=2.5)

    p.add_argument("--chunk-size", type=int, default=450)
    p.add_argument("--max-retries", type=int, default=4)
    p.add_argument("--retry-sleep-s", type=float, default=6.0)

    p.add_argument("--dpi", type=int, default=260)
    return p.parse_args()


def _launch(adql: str, *, max_retries: int, sleep_s: float) -> pd.DataFrame:
    from astroquery.gaia import Gaia

    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
    Gaia.ROW_LIMIT = -1

    last_err: Exception | None = None
    for i in range(int(max_retries)):
        try:
            job = Gaia.launch_job_async(adql, dump_to_file=False)
            return job.get_results().to_pandas()
        except Exception as e:  # noqa: BLE001
            last_err = e
            time.sleep(float(sleep_s) * (1.0 + 0.6 * i))
    raise RuntimeError(f"Gaia query failed after {max_retries} attempts") from last_err


def _load_banyan_members(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    df = df.replace({None: "", np.nan: ""})
    df["Assoc"] = df["Assoc"].astype(str).str.strip()
    df = df[(df["Assoc"] != "") & (df["Assoc"] != "------")]
    df["2MASS"] = df["2MASS"].astype(str).str.strip()
    df = df[df["2MASS"] != ""]

    # Normalize 2MASS IDs to match Gaia's original_ext_source_id (strip leading 'J').
    df["tmass_id"] = df["2MASS"].str.replace(r"^J", "", regex=True)

    df["assoc_age_myr"] = df["Assoc"].map(ASSOC_AGE_MYR).astype(float)
    return df


def _gaia_for_tmass_ids(tmass_ids: list[str], *, max_retries: int, sleep_s: float) -> pd.DataFrame:
    # Gaia stores 2MASS IDs in the best-neighbour table as original_ext_source_id.
    ids_sql = ",".join([f"'{x}'" for x in tmass_ids])
    adql = f"""
    SELECT
      x.original_ext_source_id AS tmass_id,
      x.source_id AS source_id,
      g.phot_g_mean_mag AS phot_g_mean_mag,
      g.bp_rp AS bp_rp,
      g.parallax AS parallax,
      g.parallax_over_error AS parallax_over_error,
      g.ruwe AS ruwe,
      g.teff_gspphot AS teff_gspphot
    FROM gaiadr3.tmass_psc_xsc_best_neighbour AS x
    JOIN gaiadr3.gaia_source AS g
      ON g.source_id = x.source_id
    WHERE x.original_ext_source_id IN ({ids_sql})
    """
    return _launch(adql, max_retries=max_retries, sleep_s=sleep_s)


def _design_matrix(v: np.ndarray, teff: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    teff = np.asarray(teff, dtype=float)
    lt = np.log10(np.clip(teff, 1.0, None))
    return np.column_stack([np.ones_like(v), v, lt, lt**2, v * lt])


def _fit_logsnr_model(df_etc: pd.DataFrame, y_col: str) -> np.ndarray:
    v = pd.to_numeric(df_etc["sy_vmag"], errors="coerce").to_numpy()
    teff = pd.to_numeric(df_etc["st_teff"], errors="coerce").to_numpy()
    y = pd.to_numeric(df_etc[y_col], errors="coerce").to_numpy()
    m = np.isfinite(v) & np.isfinite(teff) & np.isfinite(y) & (y > 0)
    v, teff, y = v[m], teff[m], y[m]
    X = _design_matrix(v, teff)
    logy = np.log10(y)
    p1, p99 = np.percentile(logy, [1, 99])
    logy = np.clip(logy, p1, p99)
    coef, *_ = np.linalg.lstsq(X, logy, rcond=None)
    return coef


def _predict_snr(coef: np.ndarray, v: np.ndarray, teff: np.ndarray) -> np.ndarray:
    v = np.clip(np.asarray(v, dtype=float), 0.0, 25.0)
    teff = np.clip(np.asarray(teff, dtype=float), 2500.0, 60000.0)
    X = _design_matrix(v, teff)
    logy = np.einsum("ij,j->i", X.astype(np.float64, copy=False), coef.astype(np.float64, copy=False))
    logy = np.clip(logy, -6.0, 8.0)
    return 10 ** logy


def _t_need(base_h: float, req: float, snr0: np.ndarray) -> np.ndarray:
    snr0 = np.asarray(snr0, dtype=float)
    t = float(base_h) * (float(req) / snr0) ** 2
    t[~np.isfinite(t) | (snr0 <= 0)] = np.inf
    return t


def _accumulate_star_mats(
    df: pd.DataFrame,
    *,
    dmax: float,
    vmin: float,
    vmax: float,
    dd: float,
    dv: float,
    coef_nuv: np.ndarray,
    coef_vis: np.ndarray,
    coef_nir: np.ndarray,
    req_nuv: float,
    req_vis: float,
    req_nir: float,
    base_h: float,
    t1: float,
    t5: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    dist_edges = np.arange(0.0, dmax + dd + 1e-9, dd)
    v_edges = np.arange(vmin, vmax + dv + 1e-9, dv)
    shape = (len(v_edges) - 1, len(dist_edges) - 1)

    total = np.zeros(shape, dtype=float)
    reach5 = np.zeros(shape, dtype=float)
    reach1 = np.zeros(shape, dtype=float)

    dist = pd.to_numeric(df["distance_pc"], errors="coerce").to_numpy()
    v = pd.to_numeric(df["v_est"], errors="coerce").to_numpy()
    teff = pd.to_numeric(df["teff_gspphot"], errors="coerce").to_numpy()
    m = np.isfinite(dist) & np.isfinite(v) & np.isfinite(teff)
    dist, v, teff = dist[m], v[m], teff[m]

    snr_nuv = _predict_snr(coef_nuv, v, teff)
    snr_vis = _predict_snr(coef_vis, v, teff)
    snr_nir = _predict_snr(coef_nir, v, teff)
    t_need = np.maximum.reduce(
        [
            _t_need(base_h, req_nuv, snr_nuv),
            _t_need(base_h, req_vis, snr_vis),
            _t_need(base_h, req_nir, snr_nir),
        ]
    )

    di = np.floor(dist / dd).astype(int)
    vi = np.floor((v - vmin) / dv).astype(int)
    ok = (dist >= 0) & (dist <= dmax) & (v >= vmin) & (v <= vmax) & (di >= 0) & (di < shape[1]) & (vi >= 0) & (vi < shape[0])
    di, vi, t_need = di[ok], vi[ok], t_need[ok]

    for a, b, tn in zip(vi, di, t_need, strict=False):
        total[a, b] += 1.0
        if tn <= t5:
            reach5[a, b] += 1.0
        if tn <= t1:
            reach1[a, b] += 1.0

    return dist_edges, v_edges, total, reach5, reach1


def _rgba(color_hex: str, alpha: np.ndarray) -> np.ndarray:
    import matplotlib.colors as mcolors

    rgb = np.array(mcolors.to_rgb(color_hex), dtype=float)
    rgba = np.zeros(alpha.shape + (4,), dtype=float)
    rgba[..., :3] = rgb
    rgba[..., 3] = alpha
    return rgba


def _overlay_masked(
    ax,
    *,
    mat: np.ndarray,
    color: str,
    dist_edges: np.ndarray,
    v_edges: np.ndarray,
    max_alpha: float,
    min_alpha: float,
    pctl: float = 99.0,
) -> None:
    m = mat.copy().astype(float)
    m[m <= 0] = np.nan
    if not np.isfinite(m).any():
        return
    logm = np.log10(m)
    scale = float(np.nanpercentile(logm, pctl))
    if not np.isfinite(scale) or scale <= 0:
        scale = float(np.nanmax(logm)) if np.isfinite(np.nanmax(logm)) else 1.0
    a = np.clip(logm / scale, 0.0, 1.0)
    a = min_alpha + (max_alpha - min_alpha) * a
    a = np.nan_to_num(a, nan=0.0, posinf=0.0, neginf=0.0)
    rgba = _rgba(color, a)
    ax.imshow(
        rgba,
        origin="lower",
        extent=[dist_edges[0], dist_edges[-1], v_edges[0], v_edges[-1]],
        aspect="auto",
        interpolation="nearest",
    )


def _denoise_mask(mask: np.ndarray, *, iters: int, keep_thresh: int) -> np.ndarray:
    m = mask.astype(np.uint8)
    for _ in range(int(iters)):
        p = np.pad(m, 1, mode="constant", constant_values=0)
        s = (
            p[:-2, :-2]
            + p[:-2, 1:-1]
            + p[:-2, 2:]
            + p[1:-1, :-2]
            + p[1:-1, 1:-1]
            + p[1:-1, 2:]
            + p[2:, :-2]
            + p[2:, 1:-1]
            + p[2:, 2:]
        )
        m = (s >= int(keep_thresh)).astype(np.uint8)
    return m.astype(bool)


def _load_or_build_cache(args: argparse.Namespace, members: pd.DataFrame) -> pd.DataFrame:
    cache_path = Path(args.cache_members_gaia)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    if cache_path.exists():
        df = pd.read_csv(cache_path)
        return df

    tmass_ids = sorted({x for x in members["tmass_id"].astype(str).tolist() if x})
    out: list[pd.DataFrame] = []
    for i in range(0, len(tmass_ids), int(args.chunk_size)):
        chunk = tmass_ids[i : i + int(args.chunk_size)]
        df_chunk = _gaia_for_tmass_ids(chunk, max_retries=int(args.max_retries), sleep_s=float(args.retry_sleep_s))
        out.append(df_chunk)
        print(f"Gaia crossmatch chunk {i//int(args.chunk_size)+1}: {len(df_chunk)} rows")

    df_gaia = pd.concat(out, ignore_index=True) if out else pd.DataFrame()
    if df_gaia.empty:
        raise RuntimeError("No Gaia rows returned for BANYAN members (crossmatch failed).")

    # Attach association metadata back onto Gaia rows (some 2MASS IDs may map to multiple Gaia sources; keep all).
    meta = members[["tmass_id", "Assoc", "assoc_age_myr"]].drop_duplicates()
    df = df_gaia.merge(meta, on="tmass_id", how="left")

    # Derived columns
    df["phot_g_mean_mag"] = pd.to_numeric(df["phot_g_mean_mag"], errors="coerce")
    df["bp_rp"] = pd.to_numeric(df["bp_rp"], errors="coerce")
    df["parallax"] = pd.to_numeric(df["parallax"], errors="coerce")
    df["parallax_over_error"] = pd.to_numeric(df["parallax_over_error"], errors="coerce")
    df["ruwe"] = pd.to_numeric(df["ruwe"], errors="coerce")
    df["teff_gspphot"] = pd.to_numeric(df["teff_gspphot"], errors="coerce")

    df["v_est"] = _v_from_g_bp_rp(df["phot_g_mean_mag"].to_numpy(), df["bp_rp"].to_numpy())
    df["distance_pc"] = 1000.0 / df["parallax"].to_numpy(dtype=float)

    df.to_csv(cache_path, index=False)
    print(f"Wrote cache: {cache_path}")
    return df


def _apply_quality_cuts(args: argparse.Namespace, df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    m = (
        np.isfinite(df["distance_pc"])
        & np.isfinite(df["v_est"])
        & np.isfinite(df["teff_gspphot"])
        & np.isfinite(df["ruwe"])
        & np.isfinite(df["parallax_over_error"])
        & np.isfinite(df["phot_g_mean_mag"])
        & (df["parallax"] > 0)
        & (df["distance_pc"] <= float(args.dmax_pc))
        & (df["phot_g_mean_mag"] <= float(args.gmax_guard))
        & (df["ruwe"] <= float(args.ruwe_max))
        & (df["parallax_over_error"] >= float(args.plxsnr_min))
        & (df["v_est"] >= float(args.vmin_plot))
        & (df["v_est"] <= float(args.vmax_plot))
    )
    return df[m].reset_index(drop=True)


def main() -> None:
    args = _parse_args()
    _set_mpl_cache()

    members = _load_banyan_members(Path(args.members_tsv))
    members["assoc_age_myr"] = pd.to_numeric(members["assoc_age_myr"], errors="coerce")

    # Build subsets by association age (unknown ages excluded).
    le100 = members[np.isfinite(members["assoc_age_myr"]) & (members["assoc_age_myr"] <= 100.0)]
    le10 = members[np.isfinite(members["assoc_age_myr"]) & (members["assoc_age_myr"] <= 10.0)]
    unknown_age = members[~np.isfinite(members["assoc_age_myr"])]
    older = members[np.isfinite(members["assoc_age_myr"]) & (members["assoc_age_myr"] > 100.0)]

    print(f"BANYAN members: total={len(members)}; ≤100Myr={len(le100)}; ≤10Myr={len(le10)}; unknown_age={len(unknown_age)}; >100Myr={len(older)}")
    if len(unknown_age):
        print("Unknown-age associations (excluded by default):", ", ".join(sorted(set(unknown_age["Assoc"].astype(str)))))

    # Gaia crossmatch/cache only needs the ≤100 Myr set (≤10 is a subset).
    df_cache = _load_or_build_cache(args, le100)
    df_cache["assoc_age_myr"] = pd.to_numeric(df_cache["assoc_age_myr"], errors="coerce")

    # Apply Gaia quality cuts
    df_q = _apply_quality_cuts(args, df_cache)
    df_q_le100 = df_q[np.isfinite(df_q["assoc_age_myr"]) & (df_q["assoc_age_myr"] <= 100.0)].reset_index(drop=True)
    df_q_le10 = df_q[np.isfinite(df_q["assoc_age_myr"]) & (df_q["assoc_age_myr"] <= 10.0)].reset_index(drop=True)

    # Fit proxy SNR models from ETC merged table
    df_etc = pd.read_csv(args.etc_merged)
    coef_nuv = _fit_logsnr_model(df_etc, "NUV_median_snr")
    coef_vis = _fit_logsnr_model(df_etc, "VIS_median_snr")
    coef_nir = _fit_logsnr_model(df_etc, "NIR_snr")

    import matplotlib.pyplot as plt  # noqa: E402
    from matplotlib.colors import LogNorm

    fig, axes = plt.subplots(1, 2, figsize=(11.69, 8.27), constrained_layout=False)  # A4 landscape
    fig.subplots_adjust(left=0.07, right=0.98, top=0.90, bottom=0.16, wspace=0.14)

    panels = [
        ("Young (≤100 Myr): BANYAN Σ members", df_q_le100),
        ("Very young (≤10 Myr): subset", df_q_le10),
    ]
    color = "#2ca02c"

    for ax, (title, dfp) in zip(axes, panels, strict=False):
        dist_edges, v_edges, total, reach5, reach1 = _accumulate_star_mats(
            dfp,
            dmax=float(args.dmax_pc),
            vmin=float(args.vmin_plot),
            vmax=float(args.vmax_plot),
            dd=float(args.bin_dist_pc),
            dv=float(args.bin_vmag),
            coef_nuv=coef_nuv,
            coef_vis=coef_vis,
            coef_nir=coef_nir,
            req_nuv=float(args.req_nuv),
            req_vis=float(args.req_vis),
            req_nir=float(args.req_nir),
            base_h=float(args.etc_base_hours),
            t1=float(args.t1_hours),
            t5=float(args.t5_hours),
        )

        total_plot = total.copy()
        total_plot[total_plot <= 0] = np.nan
        vmax_n = np.nanpercentile(total_plot, 99) if np.isfinite(total_plot).any() else 1.0
        bg_norm = LogNorm(vmin=1, vmax=max(3, float(vmax_n)))
        ax.pcolormesh(dist_edges, v_edges, total_plot, cmap="Greys", norm=bg_norm, shading="auto")

        _overlay_masked(
            ax,
            mat=reach5,
            color=color,
            dist_edges=dist_edges,
            v_edges=v_edges,
            max_alpha=0.25,
            min_alpha=0.02,
        )
        _overlay_masked(
            ax,
            mat=reach1,
            color=color,
            dist_edges=dist_edges,
            v_edges=v_edges,
            max_alpha=0.95,
            min_alpha=0.55,
        )

        m5 = _denoise_mask(reach5 > 0, iters=1, keep_thresh=2)
        m1 = _denoise_mask(reach1 > 0, iters=1, keep_thresh=3)
        try:
            ax.contour(
                (dist_edges[:-1] + dist_edges[1:]) / 2,
                (v_edges[:-1] + v_edges[1:]) / 2,
                m5.astype(float),
                levels=[0.5],
                colors=[color],
                linewidths=1.6,
                linestyles="--",
                alpha=0.95,
            )
            ax.contour(
                (dist_edges[:-1] + dist_edges[1:]) / 2,
                (v_edges[:-1] + v_edges[1:]) / 2,
                m1.astype(float),
                levels=[0.5],
                colors=[color],
                linewidths=2.8,
                linestyles="-",
                alpha=1.0,
            )
        except Exception:
            pass

        ax.invert_yaxis()
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("Distance (pc)")
        ax.set_ylabel("V (estimated)")
        ax.grid(alpha=0.15)
        ax.set_xlim(0, float(args.dmax_pc))
        ax.set_ylim(float(args.vmax_plot), float(args.vmin_plot))

        n_total = int(np.nansum(total))
        n_5 = int(np.nansum(reach5))
        n_1 = int(np.nansum(reach1))
        ax.text(
            0.02,
            0.98,
            f"N_total={n_total:,}; ≤{args.t5_hours:.0f}h={n_5:,}; ≤{args.t1_hours:.0f}h={n_1:,}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            bbox=dict(facecolor="white", alpha=0.75, edgecolor="none"),
        )
        ax.text(
            0.98,
            0.02,
            "solid=≤1h; dashed=≤5h",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
            color="0.25",
            bbox=dict(facecolor="white", alpha=0.65, edgecolor="none"),
        )

    fig.suptitle("Gaia DR3 young-star distance–V density with SNR-based reachability tiers (BANYAN Σ members)", fontsize=13, y=0.97)
    foot = (
        "Greyscale: BANYAN Σ member density (binned, after Gaia quality cuts). "
        "Color: reachable in all three bands (NUV+VIS+NIR) from ETC-trained SNR proxy; "
        "solid contour=≤1h, dashed contour=≤5h. V from Gaia (G,BP−RP; no extinction)."
    )
    fig.text(0.07, 0.06, foot, ha="left", va="bottom", fontsize=9, color="0.25", wrap=True)

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_prefix.with_suffix(".png"), dpi=int(args.dpi))
    fig.savefig(out_prefix.with_suffix(".pdf"))
    plt.close(fig)

    print(f"Wrote: {out_prefix.with_suffix('.png')}")


if __name__ == "__main__":
    main()


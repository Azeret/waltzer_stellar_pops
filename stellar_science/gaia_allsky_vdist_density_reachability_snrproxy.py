#!/usr/bin/env python3
"""
All-sky Gaia DR3: distance–V density by category, with *SNR-based* reachability tiers.

This script combines:
  (1) an empirical SNR proxy model fitted to the local WALTzER ETC merged table
      (`feasibility_phoenix_tdur2p5_nobs1_merged.csv`) and
  (2) Gaia DR3 binned counts in (distance, V_est, Teff) for all-sky samples.

It then estimates, per bin, whether the science-case SNR requirements in **NUV+VIS+NIR**
can be reached within **1 h** and within **5 h**, assuming **S/N ∝ √t**.

Important: this is a proposal-level population scout. The SNR proxy is trained on the
current assembled ETC sample in this repo and should be treated as approximate.
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


def _set_mpl_cache() -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(Path("stellar_science/.mplconfig").resolve()))


def _v_expr() -> str:
    # Approximate Johnson V from Gaia (G, BP-RP); no extinction correction.
    return (
        "phot_g_mean_mag + (0.01760 + 0.006860*bp_rp + 0.1732*POWER(bp_rp,2) - 0.045858*POWER(bp_rp,3))"
    )


def _parallax_cut_mas(dmax_pc: float) -> float:
    return 1000.0 / float(dmax_pc)


def _launch(adql: str) -> pd.DataFrame:
    from astroquery.gaia import Gaia

    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
    Gaia.ROW_LIMIT = -1
    job = Gaia.launch_job_async(adql, dump_to_file=False)
    return job.get_results().to_pandas()


@dataclass(frozen=True)
class Case:
    key: str
    title: str
    dmax_pc: float
    color: str
    where: str
    req_nuv: float
    req_vis: float
    req_nir: float


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--etc-merged",
        default="stellar_science/data/etc/feasibility_phoenix_tdur2p5_nobs1_merged.csv",
        help="ETC merged table used to fit the SNR proxy model.",
    )
    p.add_argument("--out-prefix", default="stellar_science/fig_gaia_allsky_vdist_density_reach_snrproxy_A4")
    p.add_argument(
        "--out-md",
        default="",
        help="Optional markdown summary output path. If empty, no .md is written.",
    )
    p.add_argument("--cache-dir", default="stellar_science/data/gaia_hist_snrproxy")

    p.add_argument("--dmax-nearby", type=float, default=500.0)
    p.add_argument("--dmax-massive", type=float, default=2000.0)

    p.add_argument("--bin-dist-pc", type=float, default=10.0)
    p.add_argument("--bin-vmag", type=float, default=0.25)
    p.add_argument("--bin-teff-k", type=float, default=200.0)
    p.add_argument("--vmin-plot", type=float, default=4.0)
    p.add_argument("--vmax-plot", type=float, default=16.0)

    # Gaia quality cuts (tunable; these control "all Gaia targets" in a robust way)
    p.add_argument("--gmax-guard", type=float, default=20.0)
    p.add_argument("--ruwe-max", type=float, default=1.4)
    p.add_argument("--plxsnr-min", type=float, default=10.0)

    # Category selection proxies (Gaia DR3 astrophysical params)
    p.add_argument("--teff-massive-min", type=float, default=15000.0)
    p.add_argument("--teff-quiet-min", type=float, default=5200.0)
    p.add_argument("--teff-quiet-max", type=float, default=6500.0)
    p.add_argument("--teff-km-max", type=float, default=5200.0)
    p.add_argument("--logg-dwarf-min", type=float, default=4.0)
    p.add_argument("--logg-pms-min", type=float, default=3.2)
    p.add_argument("--logg-pms-max", type=float, default=4.0)

    # SNR requirements (match the feasibility/pops summary defaults in this repo)
    p.add_argument("--req-massive-nuv", type=float, default=50.0)
    p.add_argument("--req-massive-vis", type=float, default=100.0)
    p.add_argument("--req-massive-nir", type=float, default=30.0)

    p.add_argument("--req-quiet-nuv", type=float, default=20.0)
    p.add_argument("--req-quiet-vis", type=float, default=200.0)
    p.add_argument("--req-quiet-nir", type=float, default=30.0)

    p.add_argument("--req-active-nuv", type=float, default=20.0)
    p.add_argument("--req-active-vis", type=float, default=100.0)
    p.add_argument("--req-active-nir", type=float, default=30.0)

    p.add_argument("--req-pms-nuv", type=float, default=5.0)
    p.add_argument("--req-pms-vis", type=float, default=300.0)
    p.add_argument("--req-pms-nir", type=float, default=30.0)

    # Time thresholds to show
    p.add_argument("--t1-hours", type=float, default=1.0)
    p.add_argument("--t5-hours", type=float, default=5.0)
    p.add_argument("--etc-base-hours", type=float, default=2.5)

    p.add_argument("--dpi", type=int, default=260)
    return p.parse_args()


def _cases(args: argparse.Namespace) -> list[Case]:
    v = _v_expr()
    gguard = float(args.gmax_guard)
    ruwe = float(args.ruwe_max)
    plxsnr = float(args.plxsnr_min)

    nearby_plx = _parallax_cut_mas(float(args.dmax_nearby))
    massive_plx = _parallax_cut_mas(float(args.dmax_massive))

    base = (
        "phot_g_mean_mag IS NOT NULL AND bp_rp IS NOT NULL AND teff_gspphot IS NOT NULL AND logg_gspphot IS NOT NULL "
        "AND parallax IS NOT NULL AND parallax > 0 "
        "AND ruwe IS NOT NULL AND parallax_over_error IS NOT NULL "
        f"AND phot_g_mean_mag <= {gguard:.3f} "
        f"AND ruwe <= {ruwe:.3f} "
        f"AND parallax_over_error >= {plxsnr:.3f} "
        f"AND ({v}) BETWEEN {args.vmin_plot:.3f} AND {args.vmax_plot:.3f} "
    )

    quiet = (
        f"{base} AND parallax >= {nearby_plx:.6f} "
        f"AND teff_gspphot BETWEEN {args.teff_quiet_min:.1f} AND {args.teff_quiet_max:.1f} "
        f"AND logg_gspphot >= {args.logg_dwarf_min:.2f}"
    )
    active = (
        f"{base} AND parallax >= {nearby_plx:.6f} "
        f"AND teff_gspphot < {args.teff_km_max:.1f} "
        f"AND logg_gspphot >= {args.logg_dwarf_min:.2f}"
    )
    pms = (
        f"{base} AND parallax >= {nearby_plx:.6f} "
        f"AND teff_gspphot < {args.teff_km_max:.1f} "
        f"AND logg_gspphot BETWEEN {args.logg_pms_min:.2f} AND {args.logg_pms_max:.2f}"
    )
    massive = (
        f"{base} AND parallax >= {massive_plx:.6f} "
        f"AND teff_gspphot >= {args.teff_massive_min:.1f}"
    )

    return [
        Case(
            key="quiet_f_g_k",
            title="Quiet FGK dwarfs (reference library)",
            dmax_pc=float(args.dmax_nearby),
            color="#1f77b4",
            where=quiet,
            req_nuv=float(args.req_quiet_nuv),
            req_vis=float(args.req_quiet_vis),
            req_nir=float(args.req_quiet_nir),
        ),
        Case(
            key="activity_km",
            title="Active K/M dwarfs (activity/flares proxy)",
            dmax_pc=float(args.dmax_nearby),
            color="#ff7f0e",
            where=active,
            req_nuv=float(args.req_active_nuv),
            req_vis=float(args.req_active_vis),
            req_nir=float(args.req_active_nir),
        ),
        Case(
            key="pms_km_proxy",
            title="PMS K/M candidates (proxy; not youth-selected)",
            dmax_pc=float(args.dmax_nearby),
            color="#2ca02c",
            where=pms,
            req_nuv=float(args.req_pms_nuv),
            req_vis=float(args.req_pms_vis),
            req_nir=float(args.req_pms_nir),
        ),
        Case(
            key="massive_hot",
            title="Massive/hot stars (winds/binaries proxy; Teff≥15000 K)",
            dmax_pc=float(args.dmax_massive),
            color="#d62728",
            where=massive,
            req_nuv=float(args.req_massive_nuv),
            req_vis=float(args.req_massive_vis),
            req_nir=float(args.req_massive_nir),
        ),
    ]


def _hist_adql(case: Case, *, dd: float, dv: float, dt: float) -> str:
    v = _v_expr()
    dist = "(1000.0/parallax)"
    dist_bin = f"FLOOR(({dist})/{dd:.6f})*{dd:.6f}"
    v_bin = f"FLOOR(({v})/{dv:.6f})*{dv:.6f}"
    teff_bin = f"FLOOR((teff_gspphot)/{dt:.6f})*{dt:.6f}"
    return f"""
    SELECT dist_bin_pc, v_bin, teff_bin_k, COUNT(*) AS n
    FROM (
      SELECT
        {dist_bin} AS dist_bin_pc,
        {v_bin} AS v_bin,
        {teff_bin} AS teff_bin_k
      FROM gaiadr3.gaia_source
      WHERE {case.where}
    ) AS sub
    GROUP BY dist_bin_pc, v_bin, teff_bin_k
    """


def _load_or_query(case: Case, cache_dir: Path, *, dd: float, dv: float, dt: float) -> pd.DataFrame:
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{case.key}_dd{dd:g}_dv{dv:g}_dt{dt:g}.csv"
    if cache_path.exists():
        return pd.read_csv(cache_path)
    df = _launch(_hist_adql(case, dd=dd, dv=dv, dt=dt))
    df.to_csv(cache_path, index=False)
    return df


def _design_matrix(v: np.ndarray, teff: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    teff = np.asarray(teff, dtype=float)
    lt = np.log10(np.clip(teff, 1.0, None))
    return np.column_stack(
        [
            np.ones_like(v),
            v,
            lt,
            lt**2,
            v * lt,
        ]
    )


def _fit_logsnr_model(df_etc: pd.DataFrame, y_col: str) -> np.ndarray:
    v = pd.to_numeric(df_etc["sy_vmag"], errors="coerce").to_numpy()
    teff = pd.to_numeric(df_etc["st_teff"], errors="coerce").to_numpy()
    y = pd.to_numeric(df_etc[y_col], errors="coerce").to_numpy()
    m = np.isfinite(v) & np.isfinite(teff) & np.isfinite(y) & (y > 0)
    v, teff, y = v[m], teff[m], y[m]
    X = _design_matrix(v, teff)
    logy = np.log10(y)

    # Simple least squares with mild clipping to reduce leverage from outliers.
    # (proposal-level proxy; keep it simple and deterministic)
    p1, p99 = np.percentile(logy, [1, 99])
    logy = np.clip(logy, p1, p99)
    coef, *_ = np.linalg.lstsq(X, logy, rcond=None)
    return coef


def _predict_snr(coef: np.ndarray, v: np.ndarray, teff: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    teff = np.asarray(teff, dtype=float)
    v = np.clip(v, 0.0, 25.0)
    teff = np.clip(teff, 2500.0, 60000.0)
    X = _design_matrix(v, teff)
    # Avoid spurious BLAS warnings by using explicit einsum (stable and fast at this size).
    logy = np.einsum("ij,j->i", X.astype(np.float64, copy=False), coef.astype(np.float64, copy=False))
    # Avoid numerical overflow in exponentiation.
    logy = np.clip(logy, -6.0, 8.0)
    return 10 ** logy


def _t_need(base_h: float, req: float, snr0: np.ndarray) -> np.ndarray:
    snr0 = np.asarray(snr0, dtype=float)
    t = float(base_h) * (float(req) / snr0) ** 2
    t[~np.isfinite(t) | (snr0 <= 0)] = np.inf
    return t


def _accumulate_mats(
    df_hist: pd.DataFrame,
    *,
    dmax: float,
    vmin: float,
    vmax: float,
    dd: float,
    dv: float,
    dt: float,
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
    # Edges for plotting
    dist_edges = np.arange(0.0, dmax + dd + 1e-9, dd)
    v_edges = np.arange(vmin, vmax + dv + 1e-9, dv)
    shape = (len(v_edges) - 1, len(dist_edges) - 1)
    total = np.zeros(shape, dtype=float)
    reach5 = np.zeros(shape, dtype=float)
    reach1 = np.zeros(shape, dtype=float)

    ddv = pd.to_numeric(df_hist["dist_bin_pc"], errors="coerce").to_numpy()
    vv = pd.to_numeric(df_hist["v_bin"], errors="coerce").to_numpy()
    tt = pd.to_numeric(df_hist["teff_bin_k"], errors="coerce").to_numpy()
    nn = pd.to_numeric(df_hist["n"], errors="coerce").to_numpy()
    good = np.isfinite(ddv) & np.isfinite(vv) & np.isfinite(tt) & np.isfinite(nn)
    ddv, vv, tt, nn = ddv[good], vv[good], tt[good], nn[good]

    # Bin centers (use +0.5 bin for predicting SNR)
    v_c = vv + 0.5 * dv
    t_c = tt + 0.5 * dt

    snr_nuv = _predict_snr(coef_nuv, v_c, t_c)
    snr_vis = _predict_snr(coef_vis, v_c, t_c)
    snr_nir = _predict_snr(coef_nir, v_c, t_c)
    t_need = np.maximum.reduce(
        [
            _t_need(base_h, req_nuv, snr_nuv),
            _t_need(base_h, req_vis, snr_vis),
            _t_need(base_h, req_nir, snr_nir),
        ]
    )

    # Map to indices
    di = np.floor(ddv / dd).astype(int)
    vi = np.floor((vv - vmin) / dv).astype(int)
    ok = (di >= 0) & (di < shape[1]) & (vi >= 0) & (vi < shape[0])
    di, vi, nn, t_need = di[ok], vi[ok], nn[ok], t_need[ok]

    for a, b, c, tn in zip(vi, di, nn, t_need, strict=False):
        total[a, b] += c
        if tn <= t5:
            reach5[a, b] += c
        if tn <= t1:
            reach1[a, b] += c

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
    """
    Overlay a colored transparency layer where mat>0, scaling alpha by log10(mat).
    """
    m = mat.copy().astype(float)
    m[m <= 0] = np.nan
    if not np.isfinite(m).any():
        return
    logm = np.log10(m)
    # Robust scaling so a few extreme bins don't wash out everything.
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
    """
    Simple 3x3 neighborhood filter to reduce fragmented contours.

    keep_thresh: minimum number of ON pixels (including self) in the 3x3 neighborhood.
    """
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


def main() -> None:
    args = _parse_args()
    _set_mpl_cache()

    # Fit proxy SNR models from ETC merged table
    df_etc = pd.read_csv(args.etc_merged)
    coef_nuv = _fit_logsnr_model(df_etc, "NUV_median_snr")
    coef_vis = _fit_logsnr_model(df_etc, "VIS_median_snr")
    coef_nir = _fit_logsnr_model(df_etc, "NIR_snr")

    cases = _cases(args)
    cache_dir = Path(args.cache_dir)

    import matplotlib.pyplot as plt  # noqa: E402
    from matplotlib.colors import LogNorm

    fig, axes = plt.subplots(2, 2, figsize=(11.69, 8.27), constrained_layout=False)  # A4 landscape
    axes = axes.ravel()
    # Give the footer its own space (avoid overlap with x-axis labels).
    fig.subplots_adjust(left=0.07, right=0.98, top=0.92, bottom=0.18, wspace=0.18, hspace=0.30)

    md_lines: list[str] = []
    md_lines.append("# Gaia DR3 all-sky distance–V density with SNR-based reachability tiers (proxy)\n\n")
    md_lines.append(f"Reachability uses an empirical SNR proxy fitted to `{args.etc_merged}` (baseline {args.etc_base_hours:.1f} h; S/N ∝ √t).\n\n")
    md_lines.append("**Massive/hot definition (proxy):** Gaia `teff_gspphot ≥ 15000 K` (proxy for O/B stars and similarly hot objects).\n\n")

    for ax, case in zip(axes, cases, strict=False):
        df_hist = _load_or_query(
            case,
            cache_dir,
            dd=float(args.bin_dist_pc),
            dv=float(args.bin_vmag),
            dt=float(args.bin_teff_k),
        )

        dist_edges, v_edges, total, reach5, reach1 = _accumulate_mats(
            df_hist,
            dmax=case.dmax_pc,
            vmin=float(args.vmin_plot),
            vmax=float(args.vmax_plot),
            dd=float(args.bin_dist_pc),
            dv=float(args.bin_vmag),
            dt=float(args.bin_teff_k),
            coef_nuv=coef_nuv,
            coef_vis=coef_vis,
            coef_nir=coef_nir,
            req_nuv=case.req_nuv,
            req_vis=case.req_vis,
            req_nir=case.req_nir,
            base_h=float(args.etc_base_hours),
            t1=float(args.t1_hours),
            t5=float(args.t5_hours),
        )

        total_plot = total.copy()
        total_plot[total_plot <= 0] = np.nan
        vmax_n = np.nanpercentile(total_plot, 99) if np.isfinite(total_plot).any() else 1.0
        bg_norm = LogNorm(vmin=1, vmax=max(3, float(vmax_n)))
        ax.pcolormesh(dist_edges, v_edges, total_plot, cmap="Greys", norm=bg_norm, shading="auto")

        # Overlay reachability: <=5h (lighter) and <=1h (darker) in the panel's category color.
        _overlay_masked(
            ax,
            mat=reach5,
            color=case.color,
            dist_edges=dist_edges,
            v_edges=v_edges,
            max_alpha=0.25,
            min_alpha=0.02,
        )
        _overlay_masked(
            ax,
            mat=reach1,
            color=case.color,
            dist_edges=dist_edges,
            v_edges=v_edges,
            max_alpha=0.95,
            min_alpha=0.55,
        )

        # Add explicit contour outlines so reachability is visible even when printed greyscale.
        # Convert reach mats to boolean masks on the same grid.
        m5 = reach5 > 0
        m1 = reach1 > 0
        if case.key == "massive_hot":
            m5 = _denoise_mask(m5, iters=2, keep_thresh=3)
            m1 = _denoise_mask(m1, iters=2, keep_thresh=4)
        else:
            m5 = _denoise_mask(m5, iters=1, keep_thresh=2)
            m1 = _denoise_mask(m1, iters=1, keep_thresh=3)
        try:
            ax.contour(
                (dist_edges[:-1] + dist_edges[1:]) / 2,
                (v_edges[:-1] + v_edges[1:]) / 2,
                m5.astype(float),
                levels=[0.5],
                colors=[case.color],
                linewidths=1.6,
                linestyles="--",
                alpha=0.95,
            )
            ax.contour(
                (dist_edges[:-1] + dist_edges[1:]) / 2,
                (v_edges[:-1] + v_edges[1:]) / 2,
                m1.astype(float),
                levels=[0.5],
                colors=[case.color],
                linewidths=2.8,
                linestyles="-",
                alpha=1.0,
            )
        except Exception:
            # If contour fails for an edge case, keep the filled overlays.
            pass

        ax.invert_yaxis()
        ax.set_title(case.title, fontsize=11)
        ax.set_xlabel("Distance (pc)")
        ax.set_ylabel("V (estimated)")
        ax.grid(alpha=0.15)
        ax.set_xlim(0, case.dmax_pc)
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
        # Minimal tier legend inside each panel.
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

        md_lines.append(f"## {case.title}\n")
        md_lines.append(f"- Distance cap: d ≤ {case.dmax_pc:.0f} pc\n")
        md_lines.append(f"- SNR requirements (NUV/VIS/NIR): {case.req_nuv:.0f}/{case.req_vis:.0f}/{case.req_nir:.0f}\n")
        md_lines.append(f"- Counts (after Gaia quality cuts): N_total={n_total:,}; reachable ≤{args.t1_hours:.0f}h={n_1:,}; ≤{args.t5_hours:.0f}h={n_5:,}\n\n")

    fig.suptitle("All-sky Gaia DR3: distance–V density with SNR-based reachability tiers (proxy)", fontsize=13, y=0.98)
    foot = (
        "Greyscale: Gaia DR3 density (binned, after quality cuts). Color: reachable in all three bands (NUV+VIS+NIR) "
        "from ETC-trained SNR proxy; solid contour=≤1h, dashed contour=≤5h. V from Gaia (G,BP−RP; no extinction)."
    )
    fig.text(0.07, 0.06, foot, ha="left", va="bottom", fontsize=9, color="0.25", wrap=True)

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_prefix.with_suffix(".png"), dpi=int(args.dpi))
    fig.savefig(out_prefix.with_suffix(".pdf"))
    plt.close(fig)

    print(f"Wrote: {out_prefix.with_suffix('.png')}")
    if str(args.out_md).strip():
        Path(args.out_md).write_text("".join(md_lines))
        print(f"Wrote: {args.out_md}")


if __name__ == "__main__":
    main()

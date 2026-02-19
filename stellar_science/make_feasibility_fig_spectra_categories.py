#!/usr/bin/env python
"""
Figure: Example 'observed' WALTzER spectra (ETC) for stellar science categories.

Creates a compact A4-ready figure with 4 categories (rows) and two
representative targets per category (columns):
  - "Bright": smallest exposure time needed to meet the SNR requirements
  - "Limit": largest exposure time (<= `--tmax-hours`) still meeting requirements

Each panel shows a single spectrum spanning NUV + VIS with band shading, plus a
single NIR photometric point:
  - NUV spectroscopy: 0.24–0.32 µm (R≈3000)
  - VIS spectroscopy: 0.425–0.797 µm (R≈3000)
The gap (0.32–0.425 µm) has no data in the current configuration.

Inputs:
  - merged feasibility CSV (includes Teff, Vmag, parallax, median SNRs)
  - Stage-1 ETC pickle (contains per-band wl/flux/variance arrays)

Notes:
  - In this ETC, NIR is modeled as band-integrated photometry.
  - ONC extinction is not applied in the example run; ONC NUV feasibility is
    optimistic for reddened members.
  - WALTzER channels are assumed simultaneous in the ETC: the integration time
    per visit is common to NUV/VIS/NIR; SNR differs by band.
  - Variable exposure-time scaling in this script assumes SNR ∝ sqrt(t) using
    the merged table's SNR values computed at `--base-hours`.
"""

from __future__ import annotations

import argparse
import os
import pickle
from pathlib import Path

import numpy as np
import pandas as pd


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--merged-csv",
        default="stellar_science/data/etc/feasibility_phoenix_tdur2p5_nobs1_merged.csv",
        help="Merged feasibility table with cluster + ETC summary columns.",
    )
    p.add_argument(
        "--pickle",
        default="stellar_science/data/etc/waltzer_snr_phoenix_tdur2p5_nobs1.pickle",
        help="Stage-1 ETC pickle output containing per-band spectra+variance.",
    )
    p.add_argument(
        "--out",
        default="stellar_science/fig_spectra_by_category_variable_t_A4.png",
        help="Output figure path (.png or .pdf).",
    )
    p.add_argument("--vis-snr", type=float, default=100.0, help="VIS median SNR cut (most categories).")
    p.add_argument("--nuv-snr", type=float, default=50.0, help="NUV median SNR cut (NUV-required categories).")
    p.add_argument(
        "--nir-snr",
        type=float,
        default=30.0,
        help="NIR SNR requirement (band-integrated photometry in the current ETC).",
    )
    p.add_argument(
        "--quiet-vis-snr",
        type=float,
        default=200.0,
        help="VIS median SNR cut for the quiet-star/reference-library category.",
    )
    p.add_argument(
        "--pms-vis-snr",
        type=float,
        default=300.0,
        help="VIS median SNR cut for the PMS (ONC) category.",
    )
    p.add_argument(
        "--pms-nuv-snr",
        type=float,
        default=5.0,
        help="NUV median SNR cut for the PMS (ONC) category (intended for rebinned/indices use).",
    )
    p.add_argument(
        "--base-hours",
        type=float,
        default=2.5,
        help="Exposure time (hours) used to compute the SNR columns in the merged CSV.",
    )
    p.add_argument(
        "--tmax-hours",
        type=float,
        default=5.0,
        help="Maximum per-target exposure time (hours) allowed for the 'limit' example selection.",
    )
    p.add_argument(
        "--dmax-pc",
        type=float,
        default=450.0,
        help="Maximum distance (pc) allowed for all non-ONC example targets (keeps examples 'nearby').",
    )
    p.add_argument(
        "--onc-dmax-pc",
        type=float,
        default=450.0,
        help="Maximum distance (pc) allowed for ONC example targets.",
    )
    p.add_argument(
        "--massive-teff-min",
        type=float,
        default=15000.0,
        help="Minimum Teff for the massive/hot category (K).",
    )
    p.add_argument(
        "--sigma-mult",
        type=float,
        default=1.0,
        help="Multiply the 1σ envelope by this factor for visibility in the plot.",
    )
    p.add_argument(
        "--shade-alpha",
        type=float,
        default=0.18,
        help="Alpha for the uncertainty shading.",
    )
    return p.parse_args()


def _set_mpl_cache() -> None:
    # Avoid font-cache warnings on systems where ~/.matplotlib is not writable.
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")


def _add_distance_columns(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    if "parallax" in d.columns:
        plx = pd.to_numeric(d["parallax"], errors="coerce")
        # Gaia parallax is in mas; crude inversion is OK for selecting near/far examples.
        d["distance_pc"] = np.where(plx > 0, 1000.0 / plx, np.nan)
    else:
        d["distance_pc"] = np.nan
    return d


def _cluster_parallax_median(df: pd.DataFrame) -> dict[str, float]:
    out: dict[str, float] = {}
    if "parallax" not in df.columns:
        return out
    for cl, g in df.groupby("cluster"):
        plx = pd.to_numeric(g["parallax"], errors="coerce")
        plx = plx[(plx > 0) & np.isfinite(plx)]
        out[str(cl)] = float(np.nanmedian(plx)) if len(plx) else float("nan")
    return out


def _safe_distance_pc(row: pd.Series, cluster_plx_med: dict[str, float]) -> float:
    """
    Return a robust-ish distance estimate in pc.
    Prefer 1/parallax if parallax is plausible for the cluster; otherwise use
    the cluster median parallax.
    """
    plx = float(row.get("parallax", float("nan")))
    cl = str(row.get("cluster", ""))
    plx_med = float(cluster_plx_med.get(cl, float("nan")))

    def inv(p: float) -> float:
        return 1000.0 / p if p and np.isfinite(p) and p > 0 else float("nan")

    # accept individual parallax if within a factor of 2 of the cluster median
    if np.isfinite(plx) and plx > 0 and np.isfinite(plx_med) and plx_med > 0:
        if 0.5 * plx_med <= plx <= 2.0 * plx_med:
            return inv(plx)
        return inv(plx_med)
    return inv(plx) if np.isfinite(plx) and plx > 0 else inv(plx_med)


def _pick_two_by_limit(
    df: pd.DataFrame,
    cat: str,
    vis_snr: float,
    nuv_snr: float,
    nir_snr: float,
    quiet_vis_snr: float,
    pms_vis_snr: float,
    pms_nuv_snr: float,
    base_hours: float,
    tmax_hours: float,
    dmax_pc: float,
    onc_dmax_pc: float,
    massive_teff_min: float,
) -> tuple[pd.Series, pd.Series]:
    """
    Pick two representative targets per category:
      - bright: requires the least exposure time to meet the SNR requirements
      - limit: requires the most exposure time while still <= tmax_hours

    We assume SNR scales ~ sqrt(time), and use the merged table's SNR columns
    (computed at base_hours) to scale exposure times.
    """
    d = df.copy()
    if "dist_pc" not in d.columns:
        d = _add_distance_columns(d)
        d["dist_pc"] = d["distance_pc"]

    # Category-specific science-grade SNR requirements (per bin / per res element):
    # We aim to demonstrate simultaneous NUV+VIS+NIR feasibility in all panels.
    if cat == "massive_hot":
        req_vis, req_nuv, req_nir = float(vis_snr), float(nuv_snr), float(nir_snr)
        sel = (d["st_teff"] >= float(massive_teff_min))
    elif cat == "quiet_f_g_k":
        req_vis, req_nuv, req_nir = float(quiet_vis_snr), 20.0, float(nir_snr)
        sel = (d["st_teff"] >= 5200) & (d["st_teff"] < 6500)
    elif cat == "activity_km":
        req_vis, req_nuv, req_nir = float(vis_snr), 20.0, float(nir_snr)
        sel = (d["st_teff"] < 5200) & (d["cluster"].isin(["M45", "NGC 2632"]))
    elif cat == "pms_onc":
        # ONC: pick a PMS K/M example and demonstrate feasibility with a pragmatic NUV requirement
        # intended for rebinned continuum/indices (extinction is not applied in this example run).
        req_vis, req_nuv, req_nir = float(pms_vis_snr), float(pms_nuv_snr), float(nir_snr)
        sel = (d["st_teff"] < 5200) & (d["cluster"] == "ONC") & (d["dist_pc"] <= float(onc_dmax_pc))
    else:
        raise ValueError(f"Unknown category: {cat}")

    sub = d[sel].copy()
    # Keep examples "nearby" unless the category explicitly targets ONC.
    if cat != "pms_onc":
        sub = sub[sub["dist_pc"] <= float(dmax_pc)].copy()

    base_hours = float(base_hours)
    tmax_hours = float(tmax_hours)

    def _t_need_series(req: float, snr: pd.Series) -> pd.Series:
        snr = pd.to_numeric(snr, errors="coerce")
        t = base_hours * (float(req) / snr) ** 2
        t = t.where((snr > 0) & np.isfinite(t), other=float("inf"))
        return t

    t_vis = _t_need_series(req_vis, sub["VIS_median_snr"])
    t_nuv = _t_need_series(req_nuv, sub["NUV_median_snr"])
    t_nir = _t_need_series(req_nir, sub["NIR_snr"])
    sub["t_need_h"] = np.maximum(np.maximum(t_vis, t_nuv), t_nir)
    # margin > 1 means already meets all requirements at base_hours
    sub["margin"] = np.sqrt(base_hours / sub["t_need_h"])

    # Bright example: smallest required time
    bright = sub.sort_values(["t_need_h", "sy_vmag"], ascending=[True, True]).iloc[0].copy()

    # Limit example: largest required time but still <= tmax_hours
    lim_sub = sub[sub["t_need_h"] <= tmax_hours].copy()
    if len(lim_sub) == 0:
        lim_sub = sub.copy()
    lim_sorted = lim_sub.sort_values(["t_need_h", "sy_vmag"], ascending=[False, False])
    limit = lim_sorted.iloc[0].copy()
    if str(limit.get("pl_name")) == str(bright.get("pl_name")) and len(lim_sorted) > 1:
        limit = lim_sorted.iloc[1].copy()

    # Attach requirements for annotation
    bright["req_nuv"] = req_nuv
    bright["req_vis"] = req_vis
    bright["req_nir"] = req_nir
    limit["req_nuv"] = req_nuv
    limit["req_vis"] = req_vis
    limit["req_nir"] = req_nir
    return bright, limit


def _load_target_spectra(pickle_path: Path, target_name: str) -> dict:
    with pickle_path.open("rb") as f:
        d = pickle.load(f)
    if target_name not in d:
        raise KeyError(f"Target {target_name!r} not found in pickle {pickle_path}")
    return d[target_name]


def _plot_nuv_vis_nir(
    ax,
    nuv: dict,
    vis: dict,
    nir: dict,
    annotate: str,
    tdur_s: float,
    ann_side: str = "left",
    sigma_mult: float = 1.0,
    shade_alpha: float = 0.18,
    show_xticklabels: bool = True,
    ann_fontsize: float = 6.8,
) -> None:
    """
    Plot NUV + VIS spectra with a *compressed* x-axis gap and add the NIR photometric band.

    We compress the empty wavelength ranges (320–425 nm and 797–900 nm) so the
    three channels fit compactly without overlaps.
    """
    # Convert to nm for x-axis labeling
    nuv_wl_nm = nuv["wl"] * 1000.0
    vis_wl_nm = vis["wl"] * 1000.0
    nir_wl_nm = float(nir["wl"][0] * 1000.0)

    tdur_s = float(max(tdur_s, 1.0))

    # Convert from photons/bin (stored for a 1 s integration in the ETC pickle)
    # to a rate density photons/s/nm so NUV/VIS (R~3000 bins) and NIR (wide photometric
    # band) are comparable on a single y-axis.
    # and NIR (wide photometric band) are comparable on a single y-axis.
    nuv_width_nm = 2.0 * np.asarray(nuv["half_widths"], dtype=float) * 1000.0
    vis_width_nm = 2.0 * np.asarray(vis["half_widths"], dtype=float) * 1000.0
    nir_width_nm = float(2.0 * np.asarray(nir["half_widths"], dtype=float)[0] * 1000.0)

    nuv_flux = np.asarray(nuv["flux"], dtype=float) / np.maximum(nuv_width_nm, 1e-12)
    vis_flux = np.asarray(vis["flux"], dtype=float) / np.maximum(vis_width_nm, 1e-12)
    nir_flux = float(nir["flux"][0]) / max(nir_width_nm, 1e-12)

    # Uncertainty of the *rate* over an exposure time t:
    #   sigma_rate = sqrt(var_1s * t) / t = sqrt(var_1s / t)
    nuv_sig = np.sqrt(np.maximum(np.asarray(nuv["variance"], dtype=float), 0.0) / tdur_s) / np.maximum(nuv_width_nm, 1e-12)
    vis_sig = np.sqrt(np.maximum(np.asarray(vis["variance"], dtype=float), 0.0) / tdur_s) / np.maximum(vis_width_nm, 1e-12)
    nir_sig = float(np.sqrt(max(float(nir["variance"][0]), 0.0) / tdur_s)) / max(nir_width_nm, 1e-12)

    # Piecewise x mapping (axis units) with compressed gaps
    nuv0, nuv1 = 240.0, 320.0
    vis0, vis1 = 425.0, 797.0
    nir0, nir1 = 900.0, 1600.0
    # Relative widths on the x-axis (arbitrary units) to make the panel compact.
    # Keep NIR *much* narrower (it's photometry) so it doesn't dominate x-range.
    W_nuv, W_vis, W_nir = 1.00, 1.60, 0.42
    # Gaps between channels (compressed empty wavelength ranges)
    # Slightly widened to prevent tick-label overlap at the boundaries.
    gap1, gap2 = 0.040, 0.040
    s_nuv = W_nuv / (nuv1 - nuv0)
    s_vis = W_vis / (vis1 - vis0)
    s_nir = W_nir / (nir1 - nir0)

    x_nuv_start = 0.0
    x_nuv_end = W_nuv
    x_vis_start = x_nuv_end + gap1
    x_vis_end = x_vis_start + W_vis
    x_nir_start = x_vis_end + gap2
    x_nir_end = x_nir_start + W_nir

    def xmap(wl_nm):
        w = np.asarray(wl_nm, dtype=float)
        x = np.full_like(w, np.nan, dtype=float)
        m = (w >= nuv0) & (w <= nuv1)
        x[m] = x_nuv_start + (w[m] - nuv0) * s_nuv
        m = (w >= vis0) & (w <= vis1)
        x[m] = x_vis_start + (w[m] - vis0) * s_vis
        m = (w >= nir0) & (w <= nir1)
        x[m] = x_nir_start + (w[m] - nir0) * s_nir
        return x

    sigma_mult = float(max(sigma_mult, 0.0))
    shade_alpha = float(min(max(shade_alpha, 0.0), 1.0))

    # Uncertainty shading is ±(sigma_mult × 1σ), but cap it for readability (low-SNR bins can explode).
    cap = 10.0

    def _fill(x, f, s, color):
        f = np.asarray(f, dtype=float)
        s = np.asarray(s, dtype=float) * sigma_mult
        f_pos = np.maximum(f, 1e-12)
        lo = np.maximum(f_pos - s, f_pos / cap)
        hi = np.minimum(f_pos + s, f_pos * cap)
        hi = np.maximum(hi, lo * 1.0001)
        ax.fill_between(x, lo, hi, color=color, alpha=shade_alpha, linewidth=0)
        ax.plot(x, f_pos, color=color, lw=0.8)
        # Return where the cap was actually used (for marking low-SNR regions)
        capped = (f_pos - s < f_pos / cap) | (f_pos + s > f_pos * cap)
        return capped

    # Band shading and labels (in mapped coordinates)
    ax.axvspan(x_nuv_start, x_nuv_end, color="tab:blue", alpha=0.06, linewidth=0)
    ax.axvspan(x_vis_start, x_vis_end, color="tab:orange", alpha=0.05, linewidth=0)
    ax.axvspan(x_nir_start, x_nir_end, color="tab:purple", alpha=0.05, linewidth=0)
    ax.axvline(x_nuv_end, color="0.65", lw=0.6)
    ax.axvline(x_vis_start, color="0.65", lw=0.6)
    ax.axvline(x_vis_end, color="0.65", lw=0.6)
    ax.axvline(x_nir_start, color="0.65", lw=0.6)
    # Place band labels using data-x + axes-y coordinates (robust to x-range changes)
    ax.text(
        x_nuv_start + 0.015 * (x_nuv_end - x_nuv_start),
        0.96,
        "NUV",
        transform=ax.get_xaxis_transform(),
        fontsize=8,
        color="tab:blue",
        ha="left",
        va="top",
    )
    ax.text(
        x_vis_start + 0.012 * (x_vis_end - x_vis_start),
        0.96,
        "VIS",
        transform=ax.get_xaxis_transform(),
        fontsize=8,
        color="tab:orange",
        ha="left",
        va="top",
    )
    # Keep the NIR label inside the axes to avoid colliding with the per-row
    # science-case headers, but pin it high so it doesn't overlap the NIR point.
    ax.text(
        x_nir_end - 0.012 * (x_nir_end - x_nir_start),
        0.96,
        "NIR phot",
        transform=ax.get_xaxis_transform(),
        fontsize=7.5,
        color="tab:purple",
        ha="right",
        va="top",
        clip_on=True,
        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85),
    )

    # Indicate compressed gaps
    ax.text((x_nuv_end + x_vis_start) / 2.0, 0.985, "gap", transform=ax.get_xaxis_transform(), fontsize=7, color="0.45", ha="center", va="top")
    ax.text((x_vis_end + x_nir_start) / 2.0, 0.985, "gap", transform=ax.get_xaxis_transform(), fontsize=7, color="0.45", ha="center", va="top")

    x_nuv = xmap(nuv_wl_nm)
    x_vis = xmap(vis_wl_nm)
    nuv_capped = _fill(x_nuv, nuv_flux, nuv_sig, "tab:blue")
    vis_capped = _fill(x_vis, vis_flux, vis_sig, "tab:orange")

    ax.set_xlim(-0.02, x_nir_end + 0.02)
    ax.set_yscale("log")
    ax.grid(alpha=0.22)

    # Y-limits based on flux-density ranges (not variances) across all channels.
    yvals = []
    for f in (nuv_flux, vis_flux):
        f = np.asarray(f, dtype=float)
        f = f[np.isfinite(f) & (f > 0)]
        if len(f):
            yvals.append(np.nanmin(f))
            yvals.append(np.nanmax(f))
    if np.isfinite(nir_flux) and nir_flux > 0:
        yvals.append(nir_flux)
    if yvals:
        y_min = float(np.nanmin(yvals))
        y_max = float(np.nanmax(yvals))
        if y_min > 0 and y_max > y_min:
            # Leave extra headroom at the bottom for the annotation box.
            ax.set_ylim(y_min * 0.15, y_max * 3.0)

    # NIR photometric point plotted "where it belongs" on the wavelength axis.
    if np.isfinite(nir_flux) and nir_flux > 0:
        xpt = float(xmap([nir_wl_nm])[0])
        ax.errorbar(
            [xpt],
            [nir_flux],
            yerr=[nir_sig] if np.isfinite(nir_sig) else None,
            fmt="o",
            ms=3.8,
            color="tab:purple",
            ecolor="tab:purple",
            elinewidth=2.0,
            capsize=6.0,
            capthick=1.5,
            zorder=6,
        )

    # Mark where the uncertainty band needed clipping (cap) so readers see low-SNR regions.
    # We draw a thin red strip along the bottom of the panel at capped wavelengths.
    y_strip = 0.02
    if np.any(nuv_capped):
        ax.scatter(x_nuv[nuv_capped], np.full(np.sum(nuv_capped), y_strip), s=4, color="tab:red", alpha=0.55,
                   transform=ax.get_xaxis_transform(), linewidths=0)
    if np.any(vis_capped):
        ax.scatter(x_vis[vis_capped], np.full(np.sum(vis_capped), y_strip), s=4, color="tab:red", alpha=0.55,
                   transform=ax.get_xaxis_transform(), linewidths=0)

    # Custom ticks labelled in true wavelength (nm)
    tick_nm = [240, 320, 425, 600, 797, 900, 1250, 1600]
    xt = [float(xmap([t])[0]) for t in tick_nm]
    ax.set_xticks(xt)
    labels = ax.set_xticklabels([str(t) for t in tick_nm], fontsize=8)
    # Prevent label overlap at compressed gaps by aligning boundary labels away
    # from each other (320↔425, 797↔900).
    for t, lab in zip(tick_nm, labels):
        if t in (320, 797):
            lab.set_ha("right")
            lab.set_fontsize(7)
        elif t in (425, 900):
            lab.set_ha("left")
            lab.set_fontsize(7)
        else:
            lab.set_ha("center")
    if not show_xticklabels:
        ax.set_xticklabels([])
        ax.tick_params(axis="x", which="both", length=0)

    # Keep the annotation inside the axes:
    # - "left": anchor on the left
    # - "right": anchor on the right and extend leftwards
    # - "center": centered (used for the one requested panel)
    if ann_side not in {"left", "right", "center"}:
        ann_side = "left"
    if ann_side == "left":
        ax_x, ax_ha = 0.02, "left"
    elif ann_side == "right":
        # Right side, but keep it out of the NIR region (NIR starts at ~0.86 in this layout).
        ax_x, ax_ha = 0.80, "right"
    else:
        ax_x, ax_ha = 0.50, "center"
    ax.text(
        ax_x,
        0.02,
        annotate,
        transform=ax.transAxes,
        ha=ax_ha,
        va="bottom",
        fontsize=float(ann_fontsize),
        clip_on=False,
        linespacing=1.0,
        bbox=dict(boxstyle="round,pad=0.16", fc="white", ec="none", alpha=0.9),
    )


def main() -> None:
    args = _parse_args()
    _set_mpl_cache()
    import matplotlib.pyplot as plt  # noqa: E402

    merged = pd.read_csv(args.merged_csv)
    required = {"cluster", "pl_name", "sy_vmag", "st_teff", "VIS_median_snr", "NUV_median_snr", "NIR_snr", "parallax"}
    missing = required - set(merged.columns)
    if missing:
        raise SystemExit(f"Missing required columns in {args.merged_csv}: {sorted(missing)}")

    categories = [
        (
            "massive_hot",
            "Massive/hot stars\n(winds, binaries)",
            f"req NUV≥{args.nuv_snr:.0f}, VIS≥{args.vis_snr:.0f}\nNIR≥{args.nir_snr:.0f}",
        ),
        (
            "quiet_f_g_k",
            "Quiet FGK stars\n(reference library)",
            f"req NUV≥20, VIS≥{args.quiet_vis_snr:.0f}\nNIR≥{args.nir_snr:.0f}",
        ),
        (
            "activity_km",
            "Active K/M dwarfs\n(activity, flares)",
            f"req NUV≥20, VIS≥{args.vis_snr:.0f}\nNIR≥{args.nir_snr:.0f}",
        ),
        (
            "pms_onc",
            "PMS K/M stars\n(ONC; no A_V here)",
            f"req NUV≥{args.pms_nuv_snr:.0f}, VIS≥{args.pms_vis_snr:.0f}\nNIR≥{args.nir_snr:.0f}",
        ),
    ]

    picks: list[tuple[pd.Series, pd.Series]] = []
    cluster_plx_med = _cluster_parallax_median(merged)
    merged = merged.copy()
    merged["dist_pc"] = merged.apply(lambda r: _safe_distance_pc(r, cluster_plx_med), axis=1)
    for key, _, _ in categories:
        picks.append(
            _pick_two_by_limit(
                merged,
                key,
                args.vis_snr,
                args.nuv_snr,
                args.nir_snr,
                args.quiet_vis_snr,
                args.pms_vis_snr,
                args.pms_nuv_snr,
                args.base_hours,
                args.tmax_hours,
                args.dmax_pc,
                args.onc_dmax_pc,
                args.massive_teff_min,
            )
        )

    # A4 landscape sizing (inches): 11.69 x 8.27
    # Use explicit margins to avoid any top-text overlap.
    fig = plt.figure(figsize=(11.69, 8.27), constrained_layout=False)

    # Layout:
    # - for each science case: [science-case header row, spectra row]
    # Keep spectra rows tall so the annotation box fits comfortably.
    height_ratios: list[float] = []
    for _ in categories:
        # Give a bit more vertical space to the science-case header row so it
        # never touches the spectra panel content.
        height_ratios.extend([0.62, 3.00])
    gs = fig.add_gridspec(
        nrows=len(height_ratios),
        ncols=2,
        height_ratios=height_ratios,
        wspace=0.015,
        hspace=0.09,
    )

    def fmt_time(hours: float) -> str:
        if not np.isfinite(hours):
            return "n/a"
        if hours < (1.0 / 60.0):
            return f"{hours*3600:.0f} s"
        if hours < 0.1:
            return f"{hours*60:.1f} min"
        if hours < 1.0:
            return f"{hours*60:.0f} min"
        return f"{hours:.1f} h"

    for i, ((cat_key, cat_label, cat_note), (near_row, far_row)) in enumerate(zip(categories, picks)):
        # Science topic header above the plots (requested).
        ax_hdr = fig.add_subplot(gs[2 * i, :])
        ax_hdr.axis("off")
        ax_hdr.text(
            0.0,
            0.92,
            cat_label.replace("\n", " "),
            fontsize=9.6,
            weight="bold",
            ha="left",
            va="top",
            clip_on=True,
        )
        ax_hdr.text(
            0.0,
            0.08,
            cat_note.replace("\n", ", "),
            fontsize=7.8,
            color="0.35",
            ha="left",
            va="bottom",
            clip_on=True,
        )

        ax_left = None
        for j, row in enumerate([near_row, far_row]):
            target = row["pl_name"]
            spec = _load_target_spectra(Path(args.pickle), target)
            dist = _safe_distance_pc(row, cluster_plx_med)
            t_h = float(row.get("t_need_h", args.base_hours))
            # Keep the figure practical: cap at tmax for display (selection already aims for <= tmax).
            t_h = min(t_h, float(args.tmax_hours))
            # Avoid "0.0h" / ultra-short times in the figure (overheads dominate in practice).
            t_h = max(t_h, 5.0 / 60.0)  # minimum 5 min shown/assumed

            def pred(snr0: float) -> float:
                snr0 = float(snr0)
                return snr0 * np.sqrt(t_h / float(args.base_hours)) if snr0 > 0 else 0.0

            nuv_pred = pred(row["NUV_median_snr"])
            vis_pred = pred(row["VIS_median_snr"])
            nir_pred = pred(row["NIR_snr"])

            req_nuv = float(row.get("req_nuv", args.nuv_snr))
            req_vis = float(row.get("req_vis", args.vis_snr))
            req_nir = float(row.get("req_nir", args.nir_snr))
            nir_frac_pct = (100.0 / nir_pred) if nir_pred > 0 else float("nan")
            annotate = (
                f"{row['cluster']}  d≈{dist:.0f} pc  V={row['sy_vmag']:.2f}  Teff={row['st_teff']:.0f}K\n"
                f"t={fmt_time(t_h)}; req NUV≥{req_nuv:.0f} VIS≥{req_vis:.0f} NIR≥{req_nir:.0f}\n"
                f"SNR@t: NUV {nuv_pred:.0f}  VIS {vis_pred:.0f}  NIR {nir_pred:.0f}  (σ/f≈{nir_frac_pct:.3g}%)"
            )
            ax = fig.add_subplot(gs[2 * i + 1, j])
            if j == 0:
                ax_left = ax
            # Annotation placement (exactly as requested):
            # - Left column: row1 left, rows2–4 right
            # - Right column: row1 left, row2 center, rows3–4 left
            ann_side_map = {
                (0, 0): "left",
                (1, 0): "right",
                (2, 0): "right",
                (3, 0): "right",
                (0, 1): "left",
                (1, 1): "center",
                (2, 1): "right",
                (3, 1): "right",
            }
            ann_side = ann_side_map.get((i, j), "left")
            _plot_nuv_vis_nir(
                ax,
                spec["nuv"],
                spec["vis"],
                spec["nir"],
                annotate,
                tdur_s=t_h * 3600.0,
                ann_side=ann_side,
                sigma_mult=float(args.sigma_mult),
                shade_alpha=float(args.shade_alpha),
                show_xticklabels=(i == len(categories) - 1),
                ann_fontsize=(6.6 if i == 0 else 6.2),
            )
            if i < len(categories) - 1:
                ax.set_xlabel("")
            else:
                ax.set_xlabel("Wavelength (nm)")
            if j == 0:
                ax.set_ylabel("Photons s$^{-1}$ nm$^{-1}$", fontsize=7.5, labelpad=1)
                ax.tick_params(axis="y", labelsize=7.5)
            else:
                ax.set_ylabel("")
                ax.tick_params(axis="y", labelsize=7.5)

    # Reserve a dedicated top band for the column headers (keeps them from
    # colliding with the first science-case header).
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.085, top=0.90)

    # Column headers (requested: slightly lower than page top).
    fig.text(0.42, 0.905, "Bright example", ha="center", va="bottom", fontsize=10)
    fig.text(0.78, 0.905, "Limit example", ha="center", va="bottom", fontsize=10)

    fig.text(
        0.08,
        0.015,
        "Red dots: wavelength bins where ±1σ shading is clipped for display.",
        ha="left",
        va="bottom",
        fontsize=8,
        color="0.35",
    )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    if out.suffix.lower() == ".png":
        fig.savefig(out.with_suffix(".pdf"))


if __name__ == "__main__":
    main()

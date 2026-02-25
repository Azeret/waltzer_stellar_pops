#!/usr/bin/env python3
"""
Figure: Example WALTzER ETC spectra for exoplanet-target host stars (from waltzer_planets.xlsx).

This mirrors the *format* of `make_feasibility_fig_spectra_categories.py`:
4 categories (rows) × 2 examples (columns: bright vs limit), with NUV+VIS spectra
and a single NIR photometric point.

Key difference vs the older stellar-science figure:
  - The ETC pickle format changed in WALTzER ETC v0.3.4+, so we compute the
    binned spectra on the fly from the stage-1 pickle using `simulate_spectrum()`.

We pick examples based on a tuned time estimate per target assuming S/N ∝ √t:
  t_need = pl_trandur * max( (req/ SNR_1transit)^2 ) across NUV/VIS/NIR
where SNR_1transit is the per-target median SNR from the ETC table (1 transit).
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def _set_mpl_cache() -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    os.environ.setdefault("XDG_CACHE_HOME", "/tmp")
    os.environ.setdefault("ASTROPY_CACHE_DIR", "/tmp/astropy-cache")
    os.environ.setdefault("ASTROPY_CONFIG_DIR", "/tmp/astropy-config")


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--merged-csv", default="data/planets/feasibility_planets_merged.csv")
    p.add_argument("--pickle", default="data/planets/waltzer_snr_planets_phoenix_nobs1.pickle")
    p.add_argument("--out", default="fig_planets_spectra_by_category_variable_t_A4.png")
    p.add_argument(
        "--flux-space",
        choices=["detected", "throughput_corrected"],
        default="detected",
        help="Plot detected counts (throughput included) or throughput-corrected flux (stare-mode output).",
    )

    p.add_argument("--tmax-hours", type=float, default=10.0, help="Cap displayed integration time (hours).")
    p.add_argument(
        "--select-max-transits",
        type=float,
        default=5.0,
        help="When selecting examples, restrict to targets needing ≤ this many transits (S/N ∝ √N) if possible.",
    )
    p.add_argument("--min-show-minutes", type=float, default=5.0, help="Floor displayed time to this (minutes).")

    # Per-category SNR requirements (defaults match the stellar-science figures in this repo).
    p.add_argument("--req-massive-nuv", type=float, default=50.0)
    p.add_argument("--req-massive-vis", type=float, default=100.0)
    p.add_argument("--req-massive-nir", type=float, default=30.0)

    p.add_argument("--req-quiet-nuv", type=float, default=20.0)
    p.add_argument("--req-quiet-vis", type=float, default=200.0)
    p.add_argument("--req-quiet-nir", type=float, default=30.0)

    p.add_argument("--req-active-nuv", type=float, default=20.0)
    p.add_argument("--req-active-vis", type=float, default=100.0)
    p.add_argument("--req-active-nir", type=float, default=30.0)

    p.add_argument("--sigma-mult", type=float, default=1.0)
    p.add_argument("--shade-alpha", type=float, default=0.18)
    return p.parse_args()


def _simulate_spectrum_from_source():
    try:
        import waltzer_etc as w  # type: ignore

        return w.simulate_spectrum
    except Exception:
        pass

    here = Path(__file__).resolve()
    for parent in [here.parent, *here.parents]:
        src = (parent / "waltzer_etc").resolve()
        if src.exists():
            sys.path.insert(0, str(src))
            import waltzer_etc as w  # type: ignore

            return w.simulate_spectrum

    raise SystemExit(
        "Could not import `waltzer_etc`. Install it (pip) or place a source checkout "
        "named `waltzer_etc/` next to this repository."
    )


def _apply_throughput(tso: dict, band: str, wl_um: np.ndarray, flux_e: np.ndarray, err_e: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Convert throughput-corrected stare-mode flux back to detected electrons by
    multiplying by the detector response used in the ETC.
    """
    det = tso[band]
    if det["det_type"] == "spectroscopy":
        det_wl = np.asarray(det["hires_wl"], dtype=float)
        det_resp = np.asarray(det["throughput"], dtype=float)
    else:
        det_wl = np.asarray([det["wl_min"], det["wl_max"]], dtype=float)
        det_resp = np.asarray(det["throughput"], dtype=float)[:2]
    resp = np.interp(np.asarray(wl_um, dtype=float), det_wl, det_resp, left=0.0, right=0.0)
    return flux_e * resp, err_e * resp


def _time_factor(
    snr_nuv: np.ndarray,
    snr_vis: np.ndarray,
    snr_nir: np.ndarray,
    req_nuv: np.ndarray,
    req_vis: np.ndarray,
    req_nir: np.ndarray,
) -> np.ndarray:
    snr_nuv = np.asarray(snr_nuv, dtype=float)
    snr_vis = np.asarray(snr_vis, dtype=float)
    snr_nir = np.asarray(snr_nir, dtype=float)
    req_nuv = np.asarray(req_nuv, dtype=float)
    req_vis = np.asarray(req_vis, dtype=float)
    req_nir = np.asarray(req_nir, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        f = np.nanmax(
            np.column_stack([(req_nuv / snr_nuv) ** 2, (req_vis / snr_vis) ** 2, (req_nir / snr_nir) ** 2]),
            axis=1,
        )
    return np.maximum(f, 1.0)


def _pick_two_by_limit(
    df: pd.DataFrame,
    *,
    cat: str,
    select_max_transits: float,
) -> tuple[pd.Series, pd.Series]:
    sub = df[df["category"] == cat].copy()
    if len(sub) == 0:
        raise SystemExit(f"No targets for category={cat!r}")

    sel = sub[np.isfinite(sub["t_factor"]) & (sub["t_factor"] <= float(select_max_transits))].copy()
    if len(sel) == 0:
        sel = sub

    bright = sel.sort_values(["t_need_h_cont", "sy_vmag"]).iloc[0]
    # Ensure "limit" is distinct when possible.
    if len(sel) == 1:
        limit = bright
    else:
        limit = sel.sort_values(["t_need_h_cont", "sy_vmag"], ascending=[False, True]).iloc[0]
        if str(limit["pl_name"]) == str(bright["pl_name"]):
            limit = sel.sort_values(["t_need_h_cont", "sy_vmag"], ascending=[False, True]).iloc[1]
    return bright, limit


def _cap_time_hours(t_h: float, *, tmax_hours: float, min_show_minutes: float) -> float:
    t_h = float(t_h)
    if not np.isfinite(t_h):
        return float("nan")
    t_h = min(t_h, float(tmax_hours))
    t_h = max(t_h, float(min_show_minutes) / 60.0)
    return t_h


def _fmt_time(hours: float) -> str:
    if not np.isfinite(hours):
        return "n/a"
    if hours < (1.0 / 60.0):
        return f"{hours*3600:.0f} s"
    if hours < 0.1:
        return f"{hours*60:.1f} min"
    if hours < 1.0:
        return f"{hours*60:.0f} min"
    return f"{hours:.2g} h"


def _plot_nuv_vis_nir(
    ax,
    *,
    nuv,
    vis,
    nir,
    tint_s: float,
    annotate: str,
    ann_side: str,
    sigma_mult: float,
    shade_alpha: float,
    show_xticklabels: bool,
    ann_fontsize: float,
) -> None:
    """
    Plot NUV+VIS spectra and a single NIR point with uncertainty shading.

    Inputs `nuv/vis/nir` are dictionaries with:
      - wl (micron), flux (e- collected), err (e- collected), half_widths (micron).
    """
    import matplotlib.pyplot as plt  # noqa: E402

    def per_nm(band: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        wl_nm = np.asarray(band["wl"], dtype=float) * 1000.0
        hw_nm = np.asarray(band["half_widths"], dtype=float) * 1000.0
        width_nm = 2.0 * np.maximum(hw_nm, 1e-12)
        flux = np.asarray(band["flux"], dtype=float) / max(tint_s, 1e-9) / width_nm
        sig = np.asarray(band["err"], dtype=float) / max(tint_s, 1e-9) / width_nm
        return wl_nm, flux, sig

    nuv_wl_nm, nuv_f, nuv_s = per_nm(nuv)
    vis_wl_nm, vis_f, vis_s = per_nm(vis)

    nir_wl_um = np.asarray(nir["wl"], dtype=float).ravel()
    nir_hw_um = np.asarray(nir["half_widths"], dtype=float).ravel()
    nir_flux_e = np.asarray(nir["flux"], dtype=float).ravel()
    nir_err_e = np.asarray(nir["err"], dtype=float).ravel()
    nir_wl_nm = float(nir_wl_um[0]) * 1000.0
    nir_hw_nm = float(nir_hw_um[0]) * 1000.0
    nir_w_nm = 2.0 * max(nir_hw_nm, 1e-12)
    nir_f = float(nir_flux_e[0]) / max(tint_s, 1e-9) / nir_w_nm
    nir_s = float(nir_err_e[0]) / max(tint_s, 1e-9) / nir_w_nm

    # Piecewise x mapping with gaps (same visual grammar as the stellar-case figure).
    nuv0, nuv1 = float(np.nanmin(nuv_wl_nm)), float(np.nanmax(nuv_wl_nm))
    vis0, vis1 = float(np.nanmin(vis_wl_nm)), float(np.nanmax(vis_wl_nm))
    nir0, nir1 = nir_wl_nm - nir_hw_nm, nir_wl_nm + nir_hw_nm

    gap1 = 0.09
    gap2 = 0.12
    W_nuv = 0.43
    W_vis = 0.45
    W_nir = 0.12

    x_nuv_end = W_nuv
    x_vis_start = x_nuv_end + gap1
    x_vis_end = x_vis_start + W_vis
    x_nir_start = x_vis_end + gap2
    x_nir_end = x_nir_start + W_nir

    s_nuv = W_nuv / max(nuv1 - nuv0, 1e-9)
    s_vis = W_vis / max(vis1 - vis0, 1e-9)
    s_nir = W_nir / max(nir1 - nir0, 1e-9)

    def xmap(wl_nm: np.ndarray) -> np.ndarray:
        wl_nm = np.asarray(wl_nm, dtype=float)
        x = np.full_like(wl_nm, np.nan, dtype=float)
        m = (wl_nm >= nuv0) & (wl_nm <= nuv1)
        x[m] = (wl_nm[m] - nuv0) * s_nuv
        m = (wl_nm >= vis0) & (wl_nm <= vis1)
        x[m] = x_vis_start + (wl_nm[m] - vis0) * s_vis
        m = (wl_nm >= nir0) & (wl_nm <= nir1)
        x[m] = x_nir_start + (wl_nm[m] - nir0) * s_nir
        return x

    # Shading with cap for readability (mark clipped bins with red dots).
    cap_factor = 10.0

    def fill(x: np.ndarray, f: np.ndarray, s: np.ndarray, color: str):
        x = np.asarray(x, dtype=float)
        f = np.asarray(f, dtype=float)
        s = np.asarray(s, dtype=float)
        lo = f - sigma_mult * s
        hi = f + sigma_mult * s
        lo_cap = np.maximum(lo, f / cap_factor)
        hi_cap = np.minimum(hi, f * cap_factor)
        clipped = (lo_cap != lo) | (hi_cap != hi)
        ax.fill_between(x, lo_cap, hi_cap, color=color, alpha=shade_alpha, linewidth=0)
        return clipped

    x_nuv = xmap(nuv_wl_nm)
    x_vis = xmap(vis_wl_nm)

    nuv_clip = fill(x_nuv, nuv_f, nuv_s, "tab:blue")
    vis_clip = fill(x_vis, vis_f, vis_s, "tab:orange")

    ax.plot(x_nuv, nuv_f, color="tab:blue", lw=1.0)
    ax.plot(x_vis, vis_f, color="tab:orange", lw=1.0)
    ax.errorbar(xmap(np.array([nir_wl_nm])), [nir_f], yerr=[nir_s], fmt="o", ms=3.6, mfc="white", mec="k",
                ecolor="0.35", elinewidth=0.8, capsize=1.6, zorder=3)

    ax.axvspan(0, x_nuv_end, color="tab:blue", alpha=0.05, linewidth=0)
    ax.axvspan(x_vis_start, x_vis_end, color="tab:orange", alpha=0.05, linewidth=0)
    ax.axvspan(x_nir_start, x_nir_end, color="tab:green", alpha=0.04, linewidth=0)
    ax.axvline(x_nuv_end, color="0.65", lw=0.6)
    ax.axvline(x_vis_start, color="0.65", lw=0.6)
    ax.axvline(x_vis_end, color="0.65", lw=0.6)
    ax.axvline(x_nir_start, color="0.65", lw=0.6)

    ax.text(0.01, 0.985, "NUV", transform=ax.transAxes, fontsize=7.2, color="tab:blue", ha="left", va="top")
    ax.text(0.44, 0.985, "VIS", transform=ax.transAxes, fontsize=7.2, color="tab:orange", ha="left", va="top")
    ax.text(0.86, 0.985, "NIR", transform=ax.transAxes, fontsize=7.2, color="tab:green", ha="left", va="top")

    # Mark clipped bins (low-SNR regions).
    y_strip = np.nanmin(np.r_[nuv_f, vis_f, [nir_f]]) * 0.02
    if np.isfinite(y_strip) and y_strip > 0:
        if np.any(nuv_clip):
            ax.scatter(x_nuv[nuv_clip], np.full(np.sum(nuv_clip), y_strip), s=4, color="tab:red", alpha=0.55, edgecolors="none", zorder=4)
        if np.any(vis_clip):
            ax.scatter(x_vis[vis_clip], np.full(np.sum(vis_clip), y_strip), s=4, color="tab:red", alpha=0.55, edgecolors="none", zorder=4)

    # Annotation
    ha = {"left": "left", "right": "right", "center": "center"}.get(ann_side, "left")
    x0 = {"left": 0.01, "center": 0.50, "right": 0.99}.get(ann_side, 0.01)
    ax.text(
        x0,
        0.02,
        annotate,
        transform=ax.transAxes,
        fontsize=ann_fontsize,
        ha=ha,
        va="bottom",
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.75", alpha=0.92),
    )

    ax.set_xlim(0, x_nir_end)
    ax.set_yscale("log")
    ax.grid(alpha=0.18)
    ax.tick_params(axis="both", labelsize=7.5)
    if not show_xticklabels:
        ax.set_xticklabels([])
        ax.tick_params(axis="x", which="both", length=0)
    else:
        # Provide a few representative wavelength ticks.
        ticks_nm = np.array([250, 300, 450, 600, 750, 1250], dtype=float)
        ax.set_xticks(xmap(ticks_nm))
        ax.set_xticklabels([f"{int(t)}" for t in ticks_nm], fontsize=7.5)


def main() -> None:
    args = _parse_args()
    _set_mpl_cache()
    import matplotlib.pyplot as plt  # noqa: E402
    import pickle  # noqa: E402

    base_dir = Path(__file__).resolve().parent

    def resolve_path(p: str) -> Path:
        pp = Path(p)
        return pp if pp.is_absolute() else (base_dir / pp)

    merged_csv = resolve_path(str(args.merged_csv))
    pickle_path = resolve_path(str(args.pickle))
    out_path = resolve_path(str(args.out))

    merged = pd.read_csv(merged_csv)
    needed = {
        "pl_name",
        "sy_vmag",
        "st_teff",
        "pl_trandur",
        "NUV_median_snr",
        "VIS_median_snr",
        "NIR_snr",
        "category",
    }
    missing = needed - set(merged.columns)
    if missing:
        raise SystemExit(f"Missing required columns in {args.merged_csv}: {sorted(missing)}")

    # Per-category requirements (same defaults as stellar-science).
    reqs = {
        "Massive/hot host (Teff≥10000K)": (args.req_massive_nuv, args.req_massive_vis, args.req_massive_nir),
        "A/F host (6500–10000K)": (args.req_quiet_nuv, args.req_quiet_vis, args.req_quiet_nir),
        "Quiet FGK host (5200–6500K)": (args.req_quiet_nuv, args.req_quiet_vis, args.req_quiet_nir),
        "K/M host (Teff<5200K)": (args.req_active_nuv, args.req_active_vis, args.req_active_nir),
    }

    merged = merged.copy()
    req_nuv = np.zeros(len(merged), dtype=float)
    req_vis = np.zeros(len(merged), dtype=float)
    req_nir = np.zeros(len(merged), dtype=float)
    for i, cat in enumerate(merged["category"].astype(str).tolist()):
        a, b, c = reqs.get(cat, reqs["Quiet FGK host (5200–6500K)"])
        req_nuv[i], req_vis[i], req_nir[i] = float(a), float(b), float(c)
    merged["req_nuv"] = req_nuv
    merged["req_vis"] = req_vis
    merged["req_nir"] = req_nir

    merged["t_factor"] = _time_factor(
        merged["NUV_median_snr"].to_numpy(float),
        merged["VIS_median_snr"].to_numpy(float),
        merged["NIR_snr"].to_numpy(float),
        merged["req_nuv"].to_numpy(float),
        merged["req_vis"].to_numpy(float),
        merged["req_nir"].to_numpy(float),
    )
    merged["t_need_h_cont"] = merged["pl_trandur"].to_numpy(float) * merged["t_factor"].to_numpy(float)

    categories = [
        ("Massive/hot host (Teff≥10000K)", "Massive / hot hosts\n(Teff ≥ 10,000 K)"),
        ("A/F host (6500–10000K)", "A/F hosts\n(6,500–10,000 K)"),
        ("Quiet FGK host (5200–6500K)", "Quiet FGK hosts\n(5,200–6,500 K)"),
        ("K/M host (Teff<5200K)", "K/M hosts\n(Teff < 5,200 K)"),
    ]

    picks: list[tuple[pd.Series, pd.Series]] = []
    for key, _ in categories:
        picks.append(_pick_two_by_limit(merged, cat=key, select_max_transits=float(args.select_max_transits)))

    # Load stage-1 pickle once.
    with open(pickle_path, "rb") as handle:
        spectra = pickle.load(handle)

    simulate_spectrum = _simulate_spectrum_from_source()

    fig = plt.figure(figsize=(11.69, 8.27), constrained_layout=False)
    height_ratios: list[float] = []
    for _ in categories:
        height_ratios.extend([0.62, 3.00])
    gs = fig.add_gridspec(
        nrows=len(height_ratios),
        ncols=2,
        height_ratios=height_ratios,
        wspace=0.015,
        hspace=0.09,
    )

    for i, ((cat_key, cat_label), (bright, limit)) in enumerate(zip(categories, picks)):
        ax_hdr = fig.add_subplot(gs[2 * i, :])
        ax_hdr.axis("off")
        ax_hdr.text(0.0, 0.92, cat_label.replace("\n", " "), fontsize=9.6, weight="bold", ha="left", va="top", clip_on=True)
        ax_hdr.text(
            0.0,
            0.08,
            "Shows star-only ETC flux spectra; NIR is photometry (single point).",
            fontsize=7.8,
            color="0.35",
            ha="left",
            va="bottom",
            clip_on=True,
        )

        for j, row in enumerate([bright, limit]):
            name = str(row["pl_name"])
            if name not in spectra:
                raise SystemExit(f"Target {name!r} not found in pickle: {pickle_path}")
            tso = spectra[name]

            t_h = _cap_time_hours(
                float(row["t_need_h_cont"]),
                tmax_hours=float(args.tmax_hours),
                min_show_minutes=float(args.min_show_minutes),
            )
            tdur_s = t_h * 3600.0
            # Compute binned star spectrum (stare mode) at the tuned time.
            bands, wl, flux, err, widths = simulate_spectrum(
                tso,
                obs_type="stare",
                obs_dur=float(t_h),
                n_obs=1,
                efficiency=float(tso["meta"]["efficiency"]),
            )
            band_to_idx = {b: k for k, b in enumerate(bands)}

            def band_dict(b: str) -> dict:
                k = band_to_idx[b]
                wl_b = np.asarray(wl[k], dtype=float)
                flux_b = np.asarray(flux[k], dtype=float)
                err_b = np.asarray(err[k], dtype=float)
                if args.flux_space == "detected":
                    flux_b, err_b = _apply_throughput(tso, b, wl_b, flux_b, err_b)
                return {"wl": wl_b, "flux": flux_b, "err": err_b, "half_widths": np.asarray(widths[k], dtype=float)}

            nuv = band_dict("nuv")
            vis = band_dict("vis")
            nir = band_dict("nir")

            def med_snr(b: dict) -> float:
                f = np.asarray(b["flux"], dtype=float)
                e = np.asarray(b["err"], dtype=float)
                m = np.isfinite(f) & np.isfinite(e) & (e > 0)
                return float(np.nanmedian(f[m] / e[m])) if np.any(m) else float("nan")

            snr_nuv = med_snr(nuv)
            snr_vis = med_snr(vis)
            snr_nir = med_snr(nir)

            annotate = (
                f"{name}\n"
                f"V={row['sy_vmag']:.2f}  Teff={row['st_teff']:.0f}K  t≈{_fmt_time(t_h)}\n"
                f"req NUV≥{row['req_nuv']:.0f} VIS≥{row['req_vis']:.0f} NIR≥{row['req_nir']:.0f}\n"
                f"SNR@t (median): NUV {snr_nuv:.0f}  VIS {snr_vis:.0f}  NIR {snr_nir:.0f}"
            )

            ax = fig.add_subplot(gs[2 * i + 1, j])
            ann_side = "left" if j == 0 else "right"
            _plot_nuv_vis_nir(
                ax,
                nuv=nuv,
                vis=vis,
                nir=nir,
                tint_s=float(tso["meta"]["efficiency"]) * tdur_s,
                annotate=annotate,
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
            else:
                ax.set_ylabel("")

    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.085, top=0.90)
    fig.text(0.42, 0.905, "Bright example", ha="center", va="bottom", fontsize=10)
    fig.text(0.78, 0.905, "Limit example", ha="center", va="bottom", fontsize=10)
    fig.text(0.08, 0.015, "Red dots: wavelength bins where ±1σ shading is clipped for display.", ha="left", va="bottom", fontsize=8, color="0.35")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=220)
    if out_path.suffix.lower() == ".png":
        fig.savefig(out_path.with_suffix(".pdf"))


if __name__ == "__main__":
    main()

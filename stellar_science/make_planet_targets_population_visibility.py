#!/usr/bin/env python3
"""
End-to-end: WALTzER ETC update + population + visibility figures for exoplanet targets.

Inputs:
  - waltzer_planets.xlsx (NASA Exoplanet Archive-style columns)

Outputs (defaults):
  - data/planets/targets_waltzer_planets.csv          (cleaned ETC input list)
  - data/planets/waltzer_snr_planets_*.csv/.pickle    (ETC stage-1 outputs)
  - data/planets/feasibility_planets_merged.csv       (targets + ETC + visibility)
  - fig_planets_population_snr_A4.(png|pdf)
  - fig_planets_visibility_2026_A4.(png|pdf)

Notes:
  - Visibility is Sun-angle only (min Sun separation), no orbit constraints.
  - ETC run uses the local ../waltzer_etc source tree (not pip-installed).
"""

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass
from datetime import date
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd


def _set_mpl_cache() -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    os.environ.setdefault("XDG_CACHE_HOME", "/tmp")
    os.environ.setdefault("ASTROPY_CACHE_DIR", "/tmp/astropy-cache")
    os.environ.setdefault("ASTROPY_CONFIG_DIR", "/tmp/astropy-config")


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--xlsx", default="waltzer_planets.xlsx", help="Input Excel file with exoplanet targets.")
    p.add_argument("--sheet", default="", help="Optional sheet name. Default: first sheet.")

    p.add_argument("--year", type=int, default=2026, help="Year used for visibility windows.")
    p.add_argument("--min-sun-sep-deg", type=float, default=110.0, help="Minimum Sun separation (deg).")

    p.add_argument("--diam-cm", type=float, default=35.0, help="Telescope diameter (cm) for ETC stage-1.")
    p.add_argument("--eff", type=float, default=0.6, help="Duty-cycle efficiency for ETC stage-1.")
    p.add_argument("--nobs", type=int, default=1, help="Number of transits/visits (ETC stage-1 stats).")
    p.add_argument(
        "--tdur-hours",
        type=float,
        default=float("nan"),
        help="If set, override all transit durations (hours). Default: use per-target `pl_trandur`.",
    )
    p.add_argument("--sed", choices=["llmodels", "bt_settl", "phoenix"], default="phoenix", help="SED library.")
    p.add_argument("--obs-mode", choices=["transit", "stare"], default="transit", help="ETC observing mode.")

    # Population summary: transits needed to reach the (stellar-science) SNR requirements.
    # Defaults match the proposal-style values used in the stellar_science figures.
    p.add_argument("--req-massive-nuv", type=float, default=50.0)
    p.add_argument("--req-massive-vis", type=float, default=100.0)
    p.add_argument("--req-massive-nir", type=float, default=30.0)

    p.add_argument("--req-quiet-nuv", type=float, default=20.0)
    p.add_argument("--req-quiet-vis", type=float, default=200.0)
    p.add_argument("--req-quiet-nir", type=float, default=30.0)

    p.add_argument("--req-active-nuv", type=float, default=20.0)
    p.add_argument("--req-active-vis", type=float, default=100.0)
    p.add_argument("--req-active-nir", type=float, default=30.0)

    # A simple K/M “PMS-like” bucket (kept for compatibility with the stellar-case logic).
    p.add_argument("--req-pms-nuv", type=float, default=5.0)
    p.add_argument("--req-pms-vis", type=float, default=300.0)
    p.add_argument("--req-pms-nir", type=float, default=30.0)

    p.add_argument("--tier1-transits", type=int, default=1, help="Reachability tier 1: max transits.")
    p.add_argument("--tier2-transits", type=int, default=5, help="Reachability tier 2: max transits.")

    p.add_argument("--out-pop", default="fig_planets_population_snr_A4.png")
    p.add_argument("--out-vis", default="fig_planets_visibility_2026_A4.png")
    p.add_argument("--out-dir", default="data/planets", help="Directory for intermediate CSV/ETC outputs.")

    return p.parse_args()


def _read_targets_xlsx(path: Path, *, sheet: str) -> pd.DataFrame:
    if not path.exists():
        raise SystemExit(f"Missing input: {path}")
    if sheet:
        df = pd.read_excel(path, sheet_name=sheet)
    else:
        df = pd.read_excel(path)
    df = df.dropna(how="all")
    # Drop common “Unnamed: …” columns.
    df = df.loc[:, [c for c in df.columns if not str(c).startswith("Unnamed:")]]
    return df


def _clean_targets(df: pd.DataFrame) -> pd.DataFrame:
    required = ["pl_name", "pl_trandur", "st_teff", "st_rad", "st_mass", "ra", "dec", "sy_vmag"]
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise SystemExit(f"Missing required columns in xlsx: {missing}")

    out = df.copy()
    for c in required:
        if c == "pl_name":
            out[c] = out[c].astype(str).str.strip()
        else:
            out[c] = pd.to_numeric(out[c], errors="coerce")

    # Keep only rows with a usable target definition.
    m = out["pl_name"].astype(str).str.len().gt(0)
    for c in required[1:]:
        m &= np.isfinite(out[c].to_numpy(dtype=float))
    out = out[m].copy()

    # De-duplicate by target name (keep first).
    out = out.drop_duplicates(subset=["pl_name"], keep="first").reset_index(drop=True)
    return out


def _waltzer_sample_from_source() -> "callable":
    """
    Import `waltzer_etc` either from an installed package, or from a nearby
    source checkout named `waltzer_etc/`.
    """
    try:
        from waltzer_etc.sample_snr import waltzer_sample  # type: ignore

        return waltzer_sample
    except Exception:
        pass

    here = Path(__file__).resolve()
    for parent in [here.parent, *here.parents]:
        src = (parent / "waltzer_etc").resolve()
        if src.exists():
            sys.path.insert(0, str(src))
            from waltzer_etc.sample_snr import waltzer_sample  # type: ignore

            return waltzer_sample

    raise SystemExit(
        "Could not import `waltzer_etc`. Install it (pip) or place a source checkout "
        "named `waltzer_etc/` next to this repository."
    )


def _waltzer_etc_version_and_sha() -> tuple[str, str]:
    """
    Best-effort version + git SHA for the sibling ../waltzer_etc repo.
    """
    version = ""
    sha = ""
    try:
        from waltzer_etc.version import __version__  # type: ignore

        version = str(__version__)
    except Exception:
        version = ""
    try:
        # If a source checkout exists nearby, grab its SHA.
        here = Path(__file__).resolve()
        repo = None
        for parent in [here.parent, *here.parents]:
            cand = (parent / "waltzer_etc").resolve()
            if cand.exists():
                repo = cand
                break
        if repo is not None:
            sha = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        sha = ""
    return version, sha


def _is_leap_year(year: int) -> bool:
    return (year % 4 == 0) and ((year % 100 != 0) or (year % 400 == 0))


@dataclass(frozen=True)
class Visibility:
    windows: str
    days_visible: int
    intervals_doy: list[tuple[int, int]]  # inclusive, 1-based day-of-year


def _visibility_from_ra_dec(
    ra_deg: float,
    dec_deg: float,
    *,
    year: int,
    min_sun_sep_deg: float,
) -> Visibility:
    from astropy.coordinates import SkyCoord, get_sun  # noqa: E402
    from astropy.time import Time  # noqa: E402
    import astropy.units as u  # noqa: E402

    n_days = 366 if _is_leap_year(year) else 365
    start = Time(f"{year}-01-01")
    times = start + np.arange(n_days) * u.day
    sun = get_sun(times)  # array SkyCoord (GCRS)

    target = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
    # Use ICRS separation consistent with ../visibility/visibility.py.
    # NOTE: `get_sun(t).icrs` is *not* the same thing here; build an ICRS
    # coordinate explicitly from the apparent GCRS RA/Dec.
    sun_icrs = SkyCoord(ra=sun.ra, dec=sun.dec, frame="icrs")
    sep = target.separation(sun_icrs).to_value(u.deg)
    vis = sep >= float(min_sun_sep_deg)

    # Convert boolean mask to intervals (inclusive, 0-based indices).
    intervals0: list[tuple[int, int]] = []
    in_seg = False
    seg_start = 0
    for i, v in enumerate(vis):
        if v and not in_seg:
            in_seg = True
            seg_start = i
        elif (not v) and in_seg:
            intervals0.append((seg_start, i - 1))
            in_seg = False
    if in_seg:
        intervals0.append((seg_start, n_days - 1))

    if not intervals0:
        return Visibility(windows="Not visible", days_visible=0, intervals_doy=[])

    # Merge wrap-around visibility (visible at year start and year end).
    if vis[0] and vis[-1] and len(intervals0) >= 2 and intervals0[0][0] == 0 and intervals0[-1][1] == n_days - 1:
        wrap = (intervals0[-1][0], intervals0[0][1])
        mid = intervals0[1:-1]
        intervals0 = [wrap, *mid]

    # Build windows string like "MMDD-MMDD, ..."
    t0 = date(year, 1, 1)
    parts = []
    for a, b in intervals0:
        da = t0.fromordinal(t0.toordinal() + a)
        db = t0.fromordinal(t0.toordinal() + b)
        parts.append(f"{da.month:02d}{da.day:02d}-{db.month:02d}{db.day:02d}")

    days_visible = int(vis.sum())
    intervals_doy = [(a + 1, b + 1) for a, b in intervals0]
    return Visibility(windows=", ".join(parts), days_visible=days_visible, intervals_doy=intervals_doy)


def _transits_needed(
    snr_nuv: np.ndarray,
    snr_vis: np.ndarray,
    snr_nir: np.ndarray,
    *,
    req_nuv: float | np.ndarray,
    req_vis: float | np.ndarray,
    req_nir: float | np.ndarray,
) -> np.ndarray:
    snr_nuv = np.asarray(snr_nuv, dtype=float)
    snr_vis = np.asarray(snr_vis, dtype=float)
    snr_nir = np.asarray(snr_nir, dtype=float)
    req_nuv = np.asarray(req_nuv, dtype=float)
    req_vis = np.asarray(req_vis, dtype=float)
    req_nir = np.asarray(req_nir, dtype=float)

    with np.errstate(divide="ignore", invalid="ignore"):
        f_nuv = (req_nuv / snr_nuv) ** 2
        f_vis = (req_vis / snr_vis) ** 2
        f_nir = (req_nir / snr_nir) ** 2
        scale = np.nanmax(np.column_stack([f_nuv, f_vis, f_nir]), axis=1)

    # If already above all thresholds, scale<1 → 1 transit.
    n = np.where(np.isfinite(scale), np.ceil(np.maximum(scale, 1.0)), np.nan)
    return n


def _assign_category_planets(df: pd.DataFrame) -> pd.Series:
    """
    Proposal-level mapping for the exoplanet-target host stars.

    This intentionally mirrors the stellar-science bucket boundaries used in this repo,
    but without cluster membership (we only have a target list).
    """
    teff = pd.to_numeric(df["st_teff"], errors="coerce")
    cat = pd.Series(index=df.index, dtype="object")
    cat.loc[teff >= 10000] = "Massive/hot host (Teff≥10000K)"
    cat.loc[(teff >= 6500) & (teff < 10000)] = "A/F host (6500–10000K)"
    cat.loc[(teff >= 5200) & (teff < 6500)] = "Quiet FGK host (5200–6500K)"
    cat.loc[teff < 5200] = "K/M host (Teff<5200K)"
    return cat.fillna("Unknown Teff")


def _requirements_by_category(args: argparse.Namespace) -> dict[str, tuple[float, float, float]]:
    return {
        "Massive/hot host (Teff≥10000K)": (
            float(args.req_massive_nuv),
            float(args.req_massive_vis),
            float(args.req_massive_nir),
        ),
        "A/F host (6500–10000K)": (
            float(args.req_quiet_nuv),
            float(args.req_quiet_vis),
            float(args.req_quiet_nir),
        ),
        "Quiet FGK host (5200–6500K)": (
            float(args.req_quiet_nuv),
            float(args.req_quiet_vis),
            float(args.req_quiet_nir),
        ),
        "K/M host (Teff<5200K)": (
            float(args.req_active_nuv),
            float(args.req_active_vis),
            float(args.req_active_nir),
        ),
        "Unknown Teff": (
            float(args.req_quiet_nuv),
            float(args.req_quiet_vis),
            float(args.req_quiet_nir),
        ),
    }


def _save_pop_figure(
    df: pd.DataFrame,
    out: Path,
    *,
    tier1_transits: int,
    tier2_transits: int,
) -> None:
    _set_mpl_cache()
    import matplotlib.pyplot as plt  # noqa: E402
    from matplotlib.colors import LogNorm  # noqa: E402

    fig = plt.figure(figsize=(11.69, 8.27), constrained_layout=True)  # A4 landscape
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.09], width_ratios=[1.0, 1.35])
    ax_bar = fig.add_subplot(gs[0, 0])
    ax_sc = fig.add_subplot(gs[0, 1])
    ax_note = fig.add_subplot(gs[1, :])
    ax_note.axis("off")

    # Per-category population counts in reachability tiers.
    display = {
        "Massive/hot host (Teff≥10000K)": "Massive / hot host\n(Teff ≥ 10,000 K)",
        "A/F host (6500–10000K)": "A/F host\n(6,500–10,000 K)",
        "Quiet FGK host (5200–6500K)": "Quiet FGK host\n(5,200–6,500 K)",
        "K/M host (Teff<5200K)": "K/M host\n(Teff < 5,200 K)",
        "Unknown Teff": "Unknown Teff",
    }
    summary = []
    for cat, g in df.groupby("category"):
        avail = int(len(g))
        t1 = int(np.sum(np.isfinite(g["n_transits_needed"]) & (g["n_transits_needed"] <= tier1_transits)))
        t2 = int(np.sum(np.isfinite(g["n_transits_needed"]) & (g["n_transits_needed"] <= tier2_transits)))
        summary.append({"category": cat, "available": avail, f"≤{tier1_transits}": t1, f"≤{tier2_transits}": t2})
    summary_df = pd.DataFrame(summary).sort_values(["available", "category"], ascending=[False, True])

    y = np.arange(len(summary_df))
    avail = summary_df["available"].to_numpy()
    t1 = summary_df[f"≤{tier1_transits}"].to_numpy()
    t2 = summary_df[f"≤{tier2_transits}"].to_numpy()

    ax_bar.barh(y, avail, color="0.90", edgecolor="0.70", label="Available")
    ax_bar.barh(y, t2, color="tab:green", alpha=0.55, label=f"Reachable (≤{tier2_transits} transits)")
    ax_bar.barh(y, t1, color="tab:blue", alpha=0.80, label=f"Reachable (≤{tier1_transits} transit)")

    # Keep count labels inside the axis.
    xmax = max(1.0, float(np.nanmax(avail))) * 1.12
    ax_bar.set_xlim(0.0, xmax)
    for yi, c in zip(y, avail, strict=False):
        x = min(float(c) + 0.5, xmax * 0.985)
        ax_bar.text(x, yi, f"{int(c)}", va="center", fontsize=9)

    ax_bar.set_yticks(y)
    ax_bar.set_yticklabels([display.get(c, c) for c in summary_df["category"].tolist()], fontsize=10)
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("Number of targets")
    ax_bar.set_title("Population: reachable with tuned transits", fontsize=12)
    ax_bar.grid(axis="x", alpha=0.25)
    ax_bar.legend(frameon=False, fontsize=9, loc="lower right")

    # Scatter: Teff vs V, colored by *time needed* (transits × transit duration).
    t_need_h = df["t_need_h"].to_numpy(dtype=float)
    finite = np.isfinite(t_need_h) & (t_need_h > 0)
    t_clip = 50.0
    t_plot = np.clip(t_need_h, 0.2, t_clip)
    sc = ax_sc.scatter(
        df.loc[finite, "st_teff"],
        df.loc[finite, "sy_vmag"],
        c=t_plot[finite],
        s=24,
        cmap="viridis",
        norm=LogNorm(vmin=0.2, vmax=t_clip),
        alpha=0.85,
        edgecolors="none",
    )
    ax_sc.invert_yaxis()
    ax_sc.set_xlabel("Teff (K)")
    ax_sc.set_ylabel("V (mag)")
    ax_sc.set_title("Targets colored by tuned time", fontsize=12)
    ax_sc.grid(alpha=0.25)
    cb = fig.colorbar(sc, ax=ax_sc, fraction=0.046, pad=0.04)
    cb.set_label("Time needed (hours; S/N ∝ √t)")
    cb.set_ticks([0.2, 0.5, 1, 2, 5, 10, 20, 50])

    fig.suptitle("WALTzER exoplanet target hosts: population + ETC reachability", fontsize=13)
    note = (
        "Per-target time is tuned from the ETC SNRs assuming S/N ∝ √t (stacking transits).  "
        f"Reachability tiers: ≤{tier1_transits} transit and ≤{tier2_transits} transits."
    )
    ax_note.text(0.5, 0.4, note, ha="center", va="center", fontsize=9, transform=ax_note.transAxes)

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=260)
    if out.suffix.lower() == ".png":
        fig.savefig(out.with_suffix(".pdf"))


def _save_visibility_figure(df: pd.DataFrame, out: Path, *, year: int) -> None:
    _set_mpl_cache()
    import matplotlib.pyplot as plt  # noqa: E402
    from matplotlib.cm import ScalarMappable  # noqa: E402
    from matplotlib.colors import Normalize  # noqa: E402

    fig = plt.figure(figsize=(11.69, 8.27), constrained_layout=True)  # A4 landscape
    gs = fig.add_gridspec(1, 2, width_ratios=[1.6, 1.0])
    ax_map = fig.add_subplot(gs[0, 0])
    ax_hist = fig.add_subplot(gs[0, 1])

    # Month boundaries for x-axis.
    t0 = date(year, 1, 1)
    month_starts = [date(year, m, 1) for m in range(1, 13)]
    month_doy = [1 + (d.toordinal() - t0.toordinal()) for d in month_starts]
    month_labels = [d.strftime("%b") for d in month_starts]

    ra_h = (df["ra_in"].to_numpy(dtype=float) / 15.0) % 24.0
    dec = df["dec_in"].to_numpy(dtype=float)
    vdays = df["vis_days"].to_numpy(dtype=float)

    norm = Normalize(vmin=np.nanmin(dec), vmax=np.nanmax(dec))
    cmap = plt.get_cmap("coolwarm")
    sm = ScalarMappable(norm=norm, cmap=cmap)

    # Visibility segments: x=day-of-year, y=RA (hours), colored by Dec.
    for i, row in df.iterrows():
        y = float(ra_h[i])
        color = cmap(norm(float(dec[i])))
        for a, b in row["vis_intervals_doy"]:
            ax_map.plot([a, b], [y, y], color=color, lw=2.0, alpha=0.85, solid_capstyle="round")

    ax_map.set_xlim(1, 366 if _is_leap_year(year) else 365)
    ax_map.set_ylim(0, 24)
    ax_map.set_yticks(np.arange(0, 25, 2))
    ax_map.set_ylabel("RA (hours)")
    ax_map.set_xlabel(f"Day of year ({year})")
    ax_map.set_title("Visibility windows (Sun-angle only)", fontsize=12)
    for x in month_doy:
        ax_map.axvline(x, color="0.88", lw=0.8, zorder=0)
    ax_map.set_xticks(month_doy)
    ax_map.set_xticklabels(month_labels, fontsize=9)

    cb = fig.colorbar(sm, ax=ax_map, fraction=0.035, pad=0.02)
    cb.set_label("Declination (deg)")

    # Histogram of visible days.
    ax_hist.hist(vdays[np.isfinite(vdays)], bins=16, color="0.25", alpha=0.85)
    ax_hist.set_xlabel("Days visible in year")
    ax_hist.set_ylabel("Number of targets")
    ax_hist.set_title("Visibility-duration distribution", fontsize=12)
    ax_hist.grid(alpha=0.25)

    fig.suptitle("WALTzER exoplanet targets: annual visibility", fontsize=13)

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=260)
    if out.suffix.lower() == ".png":
        fig.savefig(out.with_suffix(".pdf"))


def main() -> None:
    args = _parse_args()
    _set_mpl_cache()

    base_dir = Path(__file__).resolve().parent

    def resolve_path(p: str) -> Path:
        pp = Path(p)
        return pp if pp.is_absolute() else (base_dir / pp)

    xlsx_path = resolve_path(str(args.xlsx))
    out_dir = resolve_path(str(args.out_dir))
    out_pop = resolve_path(str(args.out_pop))
    out_vis = resolve_path(str(args.out_vis))

    out_dir.mkdir(parents=True, exist_ok=True)

    # 1) Read and clean targets from xlsx.
    df_xlsx = _read_targets_xlsx(xlsx_path, sheet=str(args.sheet))
    df_t = _clean_targets(df_xlsx)

    # Override transit duration if requested.
    tdur_hours = float(args.tdur_hours)
    if np.isfinite(tdur_hours):
        df_t["pl_trandur"] = float(tdur_hours)

    targets_csv = out_dir / "targets_waltzer_planets.csv"
    df_t.to_csv(targets_csv, index=False)

    # 2) Run WALTzER ETC stage-1 for these targets.
    waltzer_sample = _waltzer_sample_from_source()
    etc_out = out_dir / f"waltzer_snr_planets_{args.sed}_nobs{args.nobs}.csv"
    if np.isfinite(tdur_hours):
        etc_out = out_dir / f"waltzer_snr_planets_{args.sed}_tdur{tdur_hours:g}h_nobs{args.nobs}.csv"

    waltzer_sample(
        str(targets_csv),
        output_csv=str(etc_out),
        diameter=float(args.diam_cm),
        efficiency=float(args.eff),
        t_dur=(float(tdur_hours) if np.isfinite(tdur_hours) else None),
        n_obs=int(args.nobs),
        sed_type=str(args.sed),
        obs_mode=str(args.obs_mode),
    )

    df_etc = pd.read_csv(etc_out, comment="#")
    if "target" not in df_etc.columns:
        raise SystemExit(f"Unexpected ETC output format (missing 'target'): {etc_out}")

    # 3) Merge: xlsx metadata + ETC outputs.
    df = df_t.merge(
        df_etc,
        left_on="pl_name",
        right_on="target",
        how="left",
        validate="one_to_one",
        suffixes=("_in", "_etc"),
    )
    etc_version, etc_sha = _waltzer_etc_version_and_sha()
    df["etc_version"] = etc_version
    df["etc_git_sha"] = etc_sha

    # Category + requirements + tuned transits/time.
    df["category"] = _assign_category_planets(df)
    reqs = _requirements_by_category(args)
    req_nuv = np.zeros(len(df), dtype=float)
    req_vis = np.zeros(len(df), dtype=float)
    req_nir = np.zeros(len(df), dtype=float)
    for i, cat in enumerate(df["category"].astype(str).tolist()):
        a, b, c = reqs.get(cat, reqs["Unknown Teff"])
        req_nuv[i], req_vis[i], req_nir[i] = a, b, c
    df["req_nuv"] = req_nuv
    df["req_vis"] = req_vis
    df["req_nir"] = req_nir

    n_need = _transits_needed(
        df["NUV_median_snr"].to_numpy(dtype=float),
        df["VIS_median_snr"].to_numpy(dtype=float),
        df["NIR_snr"].to_numpy(dtype=float),
        req_nuv=req_nuv,
        req_vis=req_vis,
        req_nir=req_nir,
    )
    df["n_transits_needed"] = n_need
    df["t_need_h"] = df["pl_trandur"].to_numpy(dtype=float) * df["n_transits_needed"].to_numpy(dtype=float)

    # 4) Visibility windows for each target.
    vis = []
    for ra, dec in zip(df["ra_in"].to_numpy(dtype=float), df["dec_in"].to_numpy(dtype=float), strict=False):
        v = _visibility_from_ra_dec(
            float(ra),
            float(dec),
            year=int(args.year),
            min_sun_sep_deg=float(args.min_sun_sep_deg),
        )
        vis.append(v)
    df["vis_waltzer_year"] = [v.windows for v in vis]
    df["vis_days"] = [v.days_visible for v in vis]
    df["vis_intervals_doy"] = [v.intervals_doy for v in vis]

    merged_csv = out_dir / "feasibility_planets_merged.csv"
    df.to_csv(merged_csv, index=False)

    plan_cols = [
        "pl_name",
        "category",
        "etc_version",
        "etc_git_sha",
        "sy_vmag",
        "st_teff",
        "pl_trandur",
        "req_nuv",
        "req_vis",
        "req_nir",
        "NUV_median_snr",
        "VIS_median_snr",
        "NIR_snr",
        "n_transits_needed",
        "t_need_h",
        "vis_waltzer_year",
        "vis_days",
    ]
    plan_csv = out_dir / "observing_plan_planets.csv"
    df[plan_cols].sort_values(["category", "t_need_h", "sy_vmag"]).to_csv(plan_csv, index=False)

    # 6) Plots.
    _save_pop_figure(
        df,
        out_pop,
        tier1_transits=int(args.tier1_transits),
        tier2_transits=int(args.tier2_transits),
    )
    _save_visibility_figure(df, out_vis, year=int(args.year))

    # Console summary (kept compact).
    n = len(df)
    t1 = int(args.tier1_transits)
    t2 = int(args.tier2_transits)
    ok_t1 = np.isfinite(df["n_transits_needed"]) & (df["n_transits_needed"] <= t1)
    ok_t2 = np.isfinite(df["n_transits_needed"]) & (df["n_transits_needed"] <= t2)
    print("\nSummary:")
    print(f"- Targets: {n}")
    print(f"- Reachable within ≤{t1} transit (by category reqs): {int(ok_t1.sum())}")
    print(f"- Reachable within ≤{t2} transits (by category reqs): {int(ok_t2.sum())}")
    print(f"- Wrote: {targets_csv}")
    print(f"- Wrote: {etc_out}")
    print(f"- Wrote: {merged_csv}")
    print(f"- Wrote: {plan_csv}")
    print(f"- Wrote: {out_pop} (+.pdf if png)")
    print(f"- Wrote: {out_vis} (+.pdf if png)")


if __name__ == "__main__":
    main()

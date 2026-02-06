#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, get_sun
from astropy.time import Time
import astropy.units as u


def robust_median_mad(x: np.ndarray) -> tuple[float, float]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return float("nan"), float("nan")
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    return med, mad


def gaia_v_from_g_bp_rp(g_mag: np.ndarray, bp_rp: np.ndarray) -> np.ndarray:
    """
    Approximate Johnson V from Gaia G and (BP-RP) color.

    Uses a 3rd order polynomial relation commonly attributed to the Gaia DR2
    photometric transformations (Evans et al. 2018). This is *approximate* and
    not extinction-corrected.

    V - G ≈ 0.01760 + 0.006860*(BP-RP) + 0.1732*(BP-RP)^2 - 0.045858*(BP-RP)^3
    """
    g_mag = np.asarray(g_mag, dtype=float)
    bp_rp = np.asarray(bp_rp, dtype=float)
    x = bp_rp
    return g_mag + (0.01760 + 0.006860 * x + 0.1732 * x**2 - 0.045858 * x**3)


def visibility_windows(
    coord: SkyCoord,
    year: int = 2026,
    min_sun_sep_deg: float = 110.0,
) -> str:
    """
    Return visibility windows as 'MMDD-MMDD, ...' based on a fixed Sun-avoidance
    angle (ignores spacecraft orbit constraints).
    """
    min_sep = min_sun_sep_deg * u.deg
    start = Time(f"{year}-01-01")
    n_days = 366 if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)) else 365

    def is_visible(t: Time) -> bool:
        # Use GCRS consistently (avoids frame-conversion issues and matches
        # the logic in `visibility/visibility.py`).
        sun = get_sun(t)  # GCRS frame
        target = coord.transform_to(sun.frame)
        return target.separation(sun) >= min_sep

    stat = is_visible(start)
    ranges: list[list[str | None]] = []
    last_vis = None
    for n in range(1, n_days):
        d = start + n * u.day
        new_stat = is_visible(d)
        if new_stat == stat:
            continue
        if stat:
            d -= 1 * u.day
        mmdd = f"{d.datetime.month:02d}{d.datetime.day:02d}"
        if stat:
            if len(ranges) == 0:
                last_vis = mmdd
            else:
                ranges[-1][1] = mmdd
        else:
            ranges.append([mmdd, None])
        stat = new_stat

    if len(ranges) == 0:
        return "Not visible"
    if ranges[-1][1] is None:
        ranges[-1][1] = last_vis
    return ", ".join([f"{a}-{b}" for a, b in ranges])  # type: ignore[misc]


@dataclass(frozen=True)
class MembershipCuts:
    parallax_sigma: float = 5.0
    pm_sigma: float = 5.0
    parallax_floor_mas: float = 0.5
    pm_floor_masyr: float = 2.0
    ruwe_max: float = 1.4


def select_members(df: pd.DataFrame, cuts: MembershipCuts) -> pd.DataFrame:
    good = (
        np.isfinite(df["parallax"])
        & np.isfinite(df["pmra"])
        & np.isfinite(df["pmdec"])
        & np.isfinite(df["ruwe"])
        & (df["ruwe"] <= cuts.ruwe_max)
        & (df["parallax"] > 0)
    )
    base = df.loc[good].copy()
    if len(base) < 10:
        return base

    plx_med, plx_mad = robust_median_mad(base["parallax"].to_numpy())
    pmra_med, pmra_mad = robust_median_mad(base["pmra"].to_numpy())
    pmdec_med, pmdec_mad = robust_median_mad(base["pmdec"].to_numpy())

    plx_tol = max(cuts.parallax_floor_mas, cuts.parallax_sigma * plx_mad)
    pmra_tol = max(cuts.pm_floor_masyr, cuts.pm_sigma * pmra_mad)
    pmdec_tol = max(cuts.pm_floor_masyr, cuts.pm_sigma * pmdec_mad)

    m = (
        (np.abs(base["parallax"] - plx_med) <= plx_tol)
        & (np.abs(base["pmra"] - pmra_med) <= pmra_tol)
        & (np.abs(base["pmdec"] - pmdec_med) <= pmdec_tol)
    )
    out = base.loc[m].copy()
    out.attrs["cluster_plx_med_mas"] = plx_med
    out.attrs["cluster_pmra_med_masyr"] = pmra_med
    out.attrs["cluster_pmdec_med_masyr"] = pmdec_med
    out.attrs["cluster_plx_tol_mas"] = plx_tol
    out.attrs["cluster_pmra_tol_masyr"] = pmra_tol
    out.attrs["cluster_pmdec_tol_masyr"] = pmdec_tol
    return out


def gaia_cone_query(
    ra_deg: float,
    dec_deg: float,
    radius_deg: float,
    g_max: float,
) -> pd.DataFrame:
    from astroquery.gaia import Gaia

    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
    Gaia.ROW_LIMIT = -1

    adql = f"""
    SELECT
      designation,
      source_id,
      ra, dec,
      phot_g_mean_mag,
      phot_bp_mean_mag,
      phot_rp_mean_mag,
      bp_rp,
      teff_gspphot,
      logg_gspphot,
      parallax,
      pmra,
      pmdec,
      ruwe
    FROM gaiadr3.gaia_source
    WHERE 1=CONTAINS(
      POINT('ICRS', ra, dec),
      CIRCLE('ICRS', {ra_deg}, {dec_deg}, {radius_deg})
    )
    AND phot_g_mean_mag <= {g_max}
    AND teff_gspphot IS NOT NULL
    AND parallax IS NOT NULL
    AND pmra IS NOT NULL
    AND pmdec IS NOT NULL
    """
    job = Gaia.launch_job_async(adql, dump_to_file=False)
    table = job.get_results()
    df = table.to_pandas()
    for col in ["designation", "source_id"]:
        if col in df:
            df[col] = df[col].astype(str)
    return df


def resolve_cluster_center(name: str) -> SkyCoord:
    from astroquery.simbad import Simbad

    s = Simbad()
    # Use degree-valued coordinates (astroquery renamed 'ra(d)' -> 'ra')
    s.add_votable_fields("ra", "dec")
    r = s.query_object(name)
    if r is None or len(r) == 0:
        raise ValueError(f"Could not resolve cluster name {name!r} in SIMBAD")
    ra = float(r["ra"][0])
    dec = float(r["dec"][0])
    return SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")


def add_nearest_neighbor_metrics(
    members: pd.DataFrame,
    field: pd.DataFrame,
) -> pd.DataFrame:
    """
    Add nearest-neighbor metrics computed against *all* Gaia sources in the
    cone-search field (i.e., not only cluster members).
    """
    field_coords = SkyCoord(
        ra=field["ra"].to_numpy() * u.deg,
        dec=field["dec"].to_numpy() * u.deg,
        frame="icrs",
    )
    idx, sep2d, _ = field_coords.match_to_catalog_sky(field_coords, nthneighbor=2)
    nn_df = pd.DataFrame(
        {
            "source_id": field["source_id"].astype(str).to_numpy(),
            "nn_source_id": field["source_id"].astype(str).to_numpy()[idx],
            "nn_sep_arcsec": sep2d.to_value(u.arcsec),
            "nn_dG": field["phot_g_mean_mag"].to_numpy()[idx] - field["phot_g_mean_mag"].to_numpy(),
        }
    )
    out = members.merge(nn_df, how="left", on="source_id")
    return out


def build_targets_for_cluster(
    cluster: str,
    radius_deg: float,
    g_max: float,
    obs_hours: float,
    cuts: MembershipCuts,
    year: int,
    min_sun_sep_deg: float,
) -> pd.DataFrame:
    center = resolve_cluster_center(cluster)
    df = gaia_cone_query(center.ra.deg, center.dec.deg, radius_deg, g_max=g_max)
    if len(df) == 0:
        return df

    df["cluster"] = cluster

    members = select_members(df, cuts=cuts)
    if len(members) == 0:
        return members

    members = add_nearest_neighbor_metrics(members, field=df)

    bp_rp = members["bp_rp"].to_numpy()
    vmag = gaia_v_from_g_bp_rp(members["phot_g_mean_mag"].to_numpy(), bp_rp)
    members["sy_vmag"] = vmag

    members["pl_name"] = members["designation"]
    members["pl_trandur"] = float(obs_hours)
    members["st_teff"] = members["teff_gspphot"].astype(float)
    members["st_rad"] = np.nan
    members["st_mass"] = np.nan

    # Visibility varies smoothly across a cluster-sized field; compute once
    # for the cluster center to avoid per-target day-by-day loops.
    members[f"vis_waltzer_{year}"] = visibility_windows(
        center, year=year, min_sun_sep_deg=min_sun_sep_deg,
    )
    members["cluster_center_ra"] = center.ra.deg
    members["cluster_center_dec"] = center.dec.deg

    return members


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Build a waltz-compatible target list from Gaia DR3 cone searches around open clusters."
    )
    p.add_argument("--cluster", action="append", required=True, help="Cluster name resolvable by SIMBAD (repeatable).")
    p.add_argument("--radius-deg", type=float, default=2.0, help="Cone-search radius in degrees (default: 2.0).")
    p.add_argument("--g-max", type=float, default=12.5, help="Max Gaia G magnitude to include (default: 12.5).")
    p.add_argument("--obs-hours", type=float, default=2.5, help="pl_trandur value to write (hours). (default: 2.5)")
    p.add_argument("--year", type=int, default=2026, help="Year used for visibility windows (default: 2026).")
    p.add_argument("--min-sun-sep-deg", type=float, default=110.0, help="Sun avoidance angle (deg). (default: 110)")
    p.add_argument("--out", type=str, required=True, help="Output CSV path.")
    return p.parse_args(argv)


def main() -> None:
    args = parse_args()
    cuts = MembershipCuts()

    frames: list[pd.DataFrame] = []
    for cl in args.cluster:
        frames.append(
            build_targets_for_cluster(
                cluster=cl,
                radius_deg=args.radius_deg,
                g_max=args.g_max,
                obs_hours=args.obs_hours,
                cuts=cuts,
                year=args.year,
                min_sun_sep_deg=args.min_sun_sep_deg,
            )
        )

    if len(frames) == 0 or all(len(f) == 0 for f in frames):
        raise SystemExit("No targets returned; try increasing --radius-deg and/or --g-max.")

    df = pd.concat(frames, ignore_index=True)
    df = df.sort_values(["cluster", "sy_vmag"], ascending=[True, True]).reset_index(drop=True)

    # Put the waltz-required columns first:
    required = ["pl_name", "pl_trandur", "st_teff", "st_rad", "st_mass", "ra", "dec", "sy_vmag"]
    cols = required + [c for c in df.columns if c not in required]
    df[cols].to_csv(args.out, index=False)

    # Print a short summary:
    for cl in sorted(df["cluster"].unique()):
        n = int((df["cluster"] == cl).sum())
        print(f"{cl}: {n} targets")


if __name__ == "__main__":
    main()

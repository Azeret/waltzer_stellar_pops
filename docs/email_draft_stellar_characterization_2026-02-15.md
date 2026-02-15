Subject: WALTzER stellar characterization — cluster targets + ETC feasibility (draft summary)

Hi all,

I’ve put together a first-pass “stellar characterization” feasibility package for WALTzER, focused on nearby clusters and using the current WALTzER ETC + a simple Sun-angle visibility calculator.

What’s in place
- A shareable repo with scripts + example outputs: `Azeret/waltzer_stellar_pops`
- Reproducible target selection from Gaia DR3 (cluster cone search + basic isolation/crowding metric)
- Visibility windows (Sun-angle-only, 2026) for the initial clusters:
  - M45 (Pleiades): 2026-09-14 to 2027-01-30
  - NGC 2632 (Praesepe): 2026-11-20 to 2027-04-07
  - ONC: 2026-10-10 to 2027-02-18
- Stellar SED template download/install into the ETC expected folders:
  - PHOENIX / BT-Settl (cool stars) + LLmodels (limited grid)

ETC feasibility (illustrative, 1 visit, 2.5 h)
- Using PHOENIX photospheric templates (no extinction applied):
  - VIS is broadly feasible for a large fraction of Gaia-selected members in M45 and Praesepe at “moderate precision” thresholds.
  - NUV feasibility is much more restrictive (bright/blue/nearby objects dominate).
- Important caveat: ONC reddening is not modeled yet, so ONC NUV feasibility is optimistic; accretion-related NUV excess is also not modeled.

Where to look
- Summary + figure: `waltzer_stellar_pops/docs/feasibility_summary_phoenix_tdur2p5_nobs1_2026-02-06.md`
- Draft 2-page science theme text: `waltzer_stellar_pops/docs/science_theme_summary_stellar_characterization.md`
- Example targets + merged ETC outputs: `waltzer_stellar_pops/examples/`

Next steps (if we agree)
- Add an extinction + accretion-excess forward model for ONC PMS stars to define a robust NUV-feasible subsample.
- Add OB-star templates (e.g., TLUSTY OSTAR/BSTAR grids) for ONC massive members and rerun ETC for those.
- Expand clusters (e.g., additional nearby associations) and finalize public-survey vs consortium program split.

Best,
[Your Name]


# Feasibility summary (WALTzER ETC)

Input: `stellar_science/targets_open_clusters.csv`
ETC run: `stellar_science/waltzer_snr_phoenix_tdur2p5_nobs1.csv` (`--sed phoenix --tdur 2.5 --nobs 1 --diam 35 --eff 0.6`)
Merged table: `stellar_science/feasibility_phoenix_tdur2p5_nobs1_merged.csv`

Caution: ONC extinction/reddening is **not** applied; NUV feasibility there is likely optimistic.

## M45

Visibility (Sun-angle only, 2026): `0914-0130`
N targets: 401; V range: 5.46–12.76
- NUV depth-uncert (ppm) p16/p50/p84: 51997/364540/1174277
- VIS depth-uncert (ppm) p16/p50/p84: 3875/8821/14063
- NIR depth-uncert (ppm) p16/p50/p84: 48/135/273
- NUV median SNR counts: ≥50: 46, ≥100: 33, ≥200: 21, ≥500: 18
- VIS median SNR counts: ≥50: 401, ≥100: 339, ≥200: 163, ≥500: 44
- NIR SNR counts: ≥50: 401, ≥100: 401, ≥200: 401, ≥500: 401

Top 5 by NUV precision:

| pl_name | sy_vmag | st_teff | NUV_transit_uncert | VIS_transit_uncert | NIR_transit_uncert |
| --- | --- | --- | --- | --- | --- |
| Gaia DR3 65287458566524928 | 5.46 | 11577 | 893 | 554 | 12 |
| Gaia DR3 69812945346809600 | 5.66 | 11637 | 979 | 608 | 14 |
| Gaia DR3 66453490649681280 | 6.18 | 10090 | 1245 | 772 | 17 |
| Gaia DR3 66715105696289024 | 6.32 | 10293 | 1328 | 823 | 19 |
| Gaia DR3 66786505232432512 | 6.44 | 10878 | 1408 | 872 | 20 |

## NGC 2632

Visibility (Sun-angle only, 2026): `1120-0407`
N targets: 380; V range: 6.30–12.67
- NUV depth-uncert (ppm) p16/p50/p84: 50654/313252/935864
- VIS depth-uncert (ppm) p16/p50/p84: 3760/9247/14438
- NIR depth-uncert (ppm) p16/p50/p84: 56/159/287
- NUV median SNR counts: ≥50: 52, ≥100: 35, ≥200: 10, ≥500: 3
- VIS median SNR counts: ≥50: 380, ≥100: 309, ≥200: 127, ≥500: 43
- NIR SNR counts: ≥50: 380, ≥100: 380, ≥200: 380, ≥500: 380

Top 5 by NUV precision:

| pl_name | sy_vmag | st_teff | NUV_transit_uncert | VIS_transit_uncert | NIR_transit_uncert |
| --- | --- | --- | --- | --- | --- |
| Gaia DR3 661322060465742336 | 6.62 | 9999 | 1526 | 945 | 21 |
| Gaia DR3 664385059341602560 | 7.48 | 9705 | 2284 | 1409 | 32 |
| Gaia DR3 661290758744228224 | 7.54 | 9572 | 2349 | 1448 | 33 |
| Gaia DR3 661291136701346304 | 6.30 | 7873 | 3202 | 802 | 13 |
| Gaia DR3 661311443306610688 | 6.65 | 7290 | 3792 | 943 | 15 |

## ONC

Visibility (Sun-angle only, 2026): `1010-0218`
N targets: 202; V range: 4.76–12.87
- NUV depth-uncert (ppm) p16/p50/p84: 20517/227330/970496
- VIS depth-uncert (ppm) p16/p50/p84: 3858/8169/13996
- NIR depth-uncert (ppm) p16/p50/p84: 66/133/272
- NUV median SNR counts: ≥50: 38, ≥100: 28, ≥200: 25, ≥500: 11
- VIS median SNR counts: ≥50: 202, ≥100: 171, ≥200: 80, ≥500: 20
- NIR SNR counts: ≥50: 202, ≥100: 202, ≥200: 202, ≥500: 202

Top 5 by NUV precision:

| pl_name | sy_vmag | st_teff | NUV_transit_uncert | VIS_transit_uncert | NIR_transit_uncert |
| --- | --- | --- | --- | --- | --- |
| Gaia DR3 3017187866581304832 | 4.76 | 23082 | 353 | 407 | 11 |
| Gaia DR3 3016904707975752448 | 5.95 | 22171 | 613 | 707 | 19 |
| Gaia DR3 3017175462715713920 | 5.76 | 19808 | 674 | 643 | 16 |
| Gaia DR3 3017271326397223424 | 6.03 | 20960 | 763 | 727 | 18 |
| Gaia DR3 3023415294281663104 | 6.18 | 14487 | 966 | 774 | 18 |

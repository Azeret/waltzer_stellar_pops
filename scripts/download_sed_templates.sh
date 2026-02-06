#!/usr/bin/env bash
set -euo pipefail

# Downloads and installs the stellar SED template libraries used by the WALTzER ETC.
#
# Installs into:
#   <waltzer_etc_dir>/waltzer_etc/data/phoenix/
#   <waltzer_etc_dir>/waltzer_etc/data/bt_settl/
#
# LLMODELS is hosted on Google Drive and is not downloaded by default.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_WALTZER_ETC_DIR="$(cd "${SCRIPT_DIR}/../../waltzer_etc" 2>/dev/null && pwd || true)"

PHOENIX_URL="https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_pheonix-models_multi_v3_synphot5.tar"
BTSETTL_URL="https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_other-spectra_multi_v2_sed.tar"
LLMODELS_URL="https://drive.google.com/file/d/1pvAs8Z7RUMJrNp-JsHunZyH2vqniUnJj/view?usp=sharing"

waltzer_etc_dir="${DEFAULT_WALTZER_ETC_DIR}"

TMP_PHOENIX="/tmp/waltzer_phoenix_synphot5.tar"
TMP_BTSETTL="/tmp/waltzer_btsettl_other_spectra.tar"
TMP_LLMODELS="/tmp/waltzer_llmodels.zip"

want_phoenix=1
want_btsettl=1
want_llmodels=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --waltzer-etc-dir)
      waltzer_etc_dir="$2"
      shift 2
      ;;
    --no-phoenix) want_phoenix=0; shift ;;
    --no-bt-settl) want_btsettl=0; shift ;;
    --llmodels) want_llmodels=1; shift ;;
    -h|--help)
      cat <<'USAGE'
Usage: scripts/download_sed_templates.sh [--waltzer-etc-dir ../waltzer_etc] [--no-phoenix] [--no-bt-settl] [--llmodels]

Downloads PHOENIX and/or BT-Settl templates from STScI/MAST Reference Atlases
and installs the needed FITS files into the WALTzER ETC data directories.

Optionally downloads the WALTzER LLMODELS grid from Google Drive (requires
`gdown`, typically installed in your python environment).
USAGE
      exit 0
      ;;
    *)
      echo "Unknown arg: $1" >&2
      exit 2
      ;;
  esac
done

if [[ -z "${waltzer_etc_dir}" || ! -d "${waltzer_etc_dir}" ]]; then
  cat >&2 <<EOF
ERROR: Could not find waltzer_etc directory.

Expected it at:
  ${DEFAULT_WALTZER_ETC_DIR}

Pass it explicitly, e.g.:
  scripts/download_sed_templates.sh --waltzer-etc-dir ../waltzer_etc --llmodels
EOF
  exit 1
fi

DEST_PHOENIX="${waltzer_etc_dir}/waltzer_etc/data/phoenix"
DEST_BTSETTL="${waltzer_etc_dir}/waltzer_etc/data/bt_settl"
DEST_LLMODELS="${waltzer_etc_dir}/waltzer_etc/data/models"

mkdir -p "${DEST_PHOENIX}" "${DEST_BTSETTL}" "${DEST_LLMODELS}"

if [[ "${want_phoenix}" -eq 1 ]]; then
  if [[ -s "${TMP_PHOENIX}" ]]; then
    echo "[PHOENIX] Using cached ${TMP_PHOENIX}"
  else
    echo "[PHOENIX] Downloading ${PHOENIX_URL}"
    curl -L --fail --retry 3 --retry-delay 2 -o "${TMP_PHOENIX}" "${PHOENIX_URL}"
  fi
  echo "[PHOENIX] Extracting phoenixm00_*.fits -> ${DEST_PHOENIX}"
  # Expected paths (as documented in waltzer_etc/README.md):
  #   grp/redcat/trds/grid/phoenix/phoenixm00/phoenixm00_*.fits
  tmp_extract="$(mktemp -d /tmp/waltzer_phoenix_extract.XXXXXX)"
  tar -xf "${TMP_PHOENIX}" -C "${tmp_extract}" 'grp/redcat/trds/grid/phoenix/phoenixm00/phoenixm00_*.fits'
  if ! find "${tmp_extract}" -maxdepth 10 -type f -name "phoenixm00_*.fits" | head -n 1 | grep -q .; then
    echo "[PHOENIX] ERROR: No phoenixm00_*.fits found in archive" >&2
    rm -rf "${tmp_extract}"
    exit 1
  fi
  find "${tmp_extract}" -maxdepth 10 -type f -name "phoenixm00_*.fits" -print0 | xargs -0 -I {} mv -f {} "${DEST_PHOENIX}/"
  rm -rf "${tmp_extract}"
fi

if [[ "${want_btsettl}" -eq 1 ]]; then
  if [[ -s "${TMP_BTSETTL}" ]]; then
    echo "[BT-Settl] Using cached ${TMP_BTSETTL}"
  else
    echo "[BT-Settl] Downloading ${BTSETTL_URL}"
    curl -L --fail --retry 3 --retry-delay 2 -o "${TMP_BTSETTL}" "${BTSETTL_URL}"
  fi
  echo "[BT-Settl] Extracting phoenixm0.0_*_5.0_2011.fits -> ${DEST_BTSETTL}"
  # Expected paths (as documented in waltzer_etc/README.md):
  #   grp/redcat/trds/source/phoenixm0.0_*_5.0_2011.fits
  tmp_extract="$(mktemp -d /tmp/waltzer_btsettl_extract.XXXXXX)"
  tar -xf "${TMP_BTSETTL}" -C "${tmp_extract}" 'grp/redcat/trds/source/phoenixm0.0_*_5.0_2011.fits'
  if ! find "${tmp_extract}" -maxdepth 10 -type f -name "phoenixm0.0_*_5.0_2011.fits" | head -n 1 | grep -q .; then
    echo "[BT-Settl] ERROR: No phoenixm0.0_*_5.0_2011.fits found in archive" >&2
    rm -rf "${tmp_extract}"
    exit 1
  fi
  find "${tmp_extract}" -maxdepth 10 -type f -name "phoenixm0.0_*_5.0_2011.fits" -print0 | xargs -0 -I {} mv -f {} "${DEST_BTSETTL}/"
  rm -rf "${tmp_extract}"
fi

if [[ "${want_llmodels}" -eq 1 ]]; then
  echo "[LLMODELS] Downloading ${LLMODELS_URL}"
  if command -v gdown >/dev/null 2>&1; then
    gdown --no-cookies --fuzzy "${LLMODELS_URL}" -O "${TMP_LLMODELS}"
  elif command -v python3 >/dev/null 2>&1; then
    python3 -m gdown --no-cookies --fuzzy "${LLMODELS_URL}" -O "${TMP_LLMODELS}"
  else
    echo "[LLMODELS] ERROR: gdown not found. Install it with: python3 -m pip install gdown" >&2
    exit 1
  fi

  echo "[LLMODELS] Extracting *.flx -> ${DEST_LLMODELS}"
  unzip -o "${TMP_LLMODELS}" '*.flx' -d "${DEST_LLMODELS}" >/dev/null
  rm -rf "${DEST_LLMODELS}/__MACOSX" || true
fi

echo "Done."

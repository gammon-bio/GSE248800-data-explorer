#!/usr/bin/env bash
# 00_download.sh — fetch and unpack the GSE248800 raw 10x matrices.
#
# GSE248800_RAW.tar (~392 MB) holds 17 per-sample *_raw_feature_bc_matrix.h5
# files (UNFILTERED — all barcodes, empty droplets included). 01_process.py maps
# each GSM to its dataset/condition and applies an empty-droplet knee filter, so
# no manual per-dataset sorting is required here.
#
#   GSM7920254-259  human_pdac        (Neg/Pos x Control/WS_PDAC/C_PDAC)
#   GSM7920260-263  mouse_pancreatic  (KPP model + controls)
#   GSM7920264-270  mouse_colorectal  (C26 model + controls)
#
# Usage:  bash data-raw/00_download.sh
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
RAW="$HERE/data/raw"
mkdir -p "$RAW"
cd "$RAW"

URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248800/suppl/GSE248800_RAW.tar"
ALT="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE248800&format=file"

if [ ! -f GSE248800_RAW.tar ]; then
  echo "[00] Downloading GSE248800_RAW.tar ..."
  curl -L -C - --retry 5 --retry-delay 5 -o GSE248800_RAW.tar "$URL" \
    || curl -L -C - --retry 5 --retry-delay 5 -o GSE248800_RAW.tar "$ALT"
fi

echo "[00] Extracting ..."
tar -xf GSE248800_RAW.tar
echo "[00] Done. $(ls *_raw_feature_bc_matrix.h5 | wc -l | tr -d ' ') H5 files in $RAW"

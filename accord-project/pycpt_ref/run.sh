#!/usr/bin/env bash
# Run the pyCPT reference for one case inside the pycpt image from accord-research.
# Usage: pycpt_ref/run.sh <case name> [/path/to/accord-research/pycpt to (re)build the image]
set -euo pipefail
CASE="${1:?case name under pycpt_ref/cases}"
HERE="$(cd "$(dirname "$0")" && pwd)"
if [ -n "${2:-}" ]; then
  docker build --platform linux/amd64 --build-arg PYCPT_VERSION=2.10.4 -t pycpt:2.10.4 "$2"
fi
docker run --rm --platform linux/amd64 -e MPLBACKEND=Agg \
  -v "$HERE/cases/$CASE:/case" -v "$HERE/run_reference.py:/run_reference.py:ro" \
  pycpt:2.10.4 python /run_reference.py

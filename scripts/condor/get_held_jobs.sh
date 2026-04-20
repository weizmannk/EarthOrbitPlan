#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# bash get_held_jobs.sh uvex_cmds.sh
#
# Dump the executable path of every held HTCondor job into a file.
# Useful for diagnosing or re-running jobs that entered a held state.
#
# Usage:
#   bash get_held_jobs.sh                  # writes to uvex_cmds.sh (default)
#   bash get_held_jobs.sh my_output.sh     # custom output file
#
# Example output (one executable path per line):
#   /path/to/pipeline/run_injection_001.sh
#   /path/to/pipeline/run_injection_002.sh
#
# Dependencies: HTCondor (condor_q), grep, sed
# -----------------------------------------------------------------------------

set -euo pipefail

OUTPUT="${1:-uvex_cmds.sh}"

condor_q -hold -long \
    | grep  "^Cmd" \
    | sed   's/Cmd = "//; s/"$//' \
    > "$OUTPUT"

echo "Held job commands written to: $OUTPUT ($(wc -l < "$OUTPUT") jobs)"

#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# List the held HTCondor jobs and dump their wrapper paths to a file.
#
# This is a diagnostic aid. It is NOT the way to recover held jobs: the file it
# writes is a plain list of wrapper paths, so running it would execute every
# m4opt solve locally on the submit node rather than through Condor. To recover
# jobs held for exceeding their memory limit, edit and release them in place:
#
#     bash release_held_jobs.sh --apply
#
# Usage:
#   bash get_held_jobs.sh                   # writes to held_cmds.txt
#   bash get_held_jobs.sh my_output.txt     # custom output file
#   bash get_held_jobs.sh my_output.txt 24177   # ... only for cluster 24177
#
# Dependencies: HTCondor (condor_q)
# -----------------------------------------------------------------------------

set -euo pipefail

OUTPUT="${1:-held_cmds.txt}"
CLUSTER="${2:-}"

CONSTRAINT="JobStatus == 5"
[ -n "$CLUSTER" ] && CONSTRAINT="$CONSTRAINT && ClusterId == $CLUSTER"

# -af reads the attribute directly, so there is no ClassAd text to parse, and
# an empty queue yields an empty result instead of tripping pipefail on grep.
condor_q -constraint "$CONSTRAINT" -af Cmd > "$OUTPUT"

N=$(wc -l < "$OUTPUT" | tr -d ' ')
echo "Held jobs: $N  (wrapper paths written to $OUTPUT)"

if [ "$N" -eq 0 ]; then
    exit 0
fi

echo
echo "Why they are held:"
condor_q -constraint "$CONSTRAINT" -af HoldReason | sort | uniq -c | sort -rn

echo
echo "Memory actually used before being killed (MB):"
condor_q -constraint "$CONSTRAINT" -af MemoryUsage | sort -n | uniq -c

echo
echo "To recover them without resubmitting:  bash ${0%/*}/release_held_jobs.sh --apply"

#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Repair HTCondor jobs held for exceeding their memory limit, in place.
#
# Jobs killed with "gone over cgroup memory limit" stay in the queue as held.
# They do not need to be resubmitted: RequestMemory can be edited and the job
# released, which reruns it with the new limit.
#
# By default this installs the same adaptive request that new submissions use
# (see earthorbitplan/backend/condor.py): the job re-requests GROWTH times the
# memory it actually used before being killed, and PeriodicRelease is set so
# that any further memory hold is retried automatically without you watching
# the queue. Pass an explicit value in MB to force a flat request instead.
#
# Usage:
#   bash release_held_jobs.sh                # dry run: show what would change
#   bash release_held_jobs.sh --apply        # adaptive: 2x measured usage
#   bash release_held_jobs.sh --apply 24177  # ... only for cluster 24177
#   bash release_held_jobs.sh 120000         # flat 120000 MB (no escalation)
#
# Escalation is capped at CAP_MB, which must stay below the RAM of the largest
# slot in the pool:  condor_status -af Memory | sort -n | tail -1
#
# Dependencies: HTCondor (condor_q, condor_qedit, condor_release)
# -----------------------------------------------------------------------------

set -euo pipefail

GROWTH="${GROWTH:-2.0}"
MAX_RETRIES="${MAX_RETRIES:-3}"
FLOOR_MB="${FLOOR_MB:-40000}"
CAP_MB="${CAP_MB:-200000}"

MODE=""
MEM_MB=""
CLUSTER=""

case "${1:-}" in
    --apply) MODE="adaptive"; CLUSTER="${2:-}" ;;
    "")      MODE="dryrun" ;;
    *)       MODE="flat"; MEM_MB="$1"; CLUSTER="${2:-}" ;;
esac

CONSTRAINT="JobStatus == 5 && Owner == \"$(whoami)\""
[ -n "$CLUSTER" ] && CONSTRAINT="$CONSTRAINT && ClusterId == $CLUSTER"

N=$(condor_q -hold -constraint "$CONSTRAINT" -af ClusterId | wc -l | tr -d ' ')
echo "Held jobs matching: $N"

if [ "$N" -eq 0 ]; then
    exit 0
fi

echo
echo "Current requests, measured usage, and hold reasons:"
condor_q -hold -constraint "$CONSTRAINT" \
    -af ClusterId ProcId RequestMemory MemoryUsage HoldReason \
    | sort -u -k5 | head -5

if [ "$MODE" = "dryrun" ]; then
    echo
    echo "Dry run. Re-run to apply:"
    echo "    bash $0 --apply     # adaptive: ${GROWTH}x measured usage, auto-retry"
    echo "    bash $0 120000      # flat 120000 MB"
    exit 0
fi

if [ "$MODE" = "adaptive" ]; then
    EXPR="ifthenelse(MemoryUsage =!= undefined, min({${CAP_MB}, max({${FLOOR_MB}, MemoryUsage * ${GROWTH}})}), ${FLOOR_MB})"
    echo
    echo "Setting adaptive RequestMemory on $N job(s):"
    echo "    $EXPR"
    condor_qedit -constraint "$CONSTRAINT" RequestMemory "$EXPR"
    # Let any future memory hold retry itself instead of sitting in the queue.
    condor_qedit -constraint "$CONSTRAINT" PeriodicRelease \
        "(HoldReasonCode == 34 || HoldReasonCode == 104) && (NumJobStarts <= ${MAX_RETRIES})" \
        || echo "warning: could not set PeriodicRelease; jobs will need releasing by hand"
else
    echo
    echo "Setting flat RequestMemory = $MEM_MB MB on $N job(s)."
    condor_qedit -constraint "$CONSTRAINT" RequestMemory "$MEM_MB"
fi

echo "Releasing..."
condor_release -constraint "$CONSTRAINT"
echo "Done. Watch them with: condor_q -nobatch"

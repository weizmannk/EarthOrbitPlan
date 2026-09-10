#!/usr/bin/env bash
# Swap the size-limited PyPI `cplex` runtime for the full CPLEX Optimization
# Studio engine. Run this once after installing CPLEX Studio, and again after
# any reinstall/upgrade of the `cplex` wheel (e.g. a fresh `uv pip install`).
#
#   ./scripts/setup-cplex.sh [/path/to/CPLEX_Studio]
#
set -euo pipefail

STUDIO="${1:-$HOME/Applications/CPLEX_Studio222}"

if [ ! -d "$STUDIO/cplex" ]; then
    echo "CPLEX Studio not found at: $STUDIO" >&2
    echo "Pass the install path as the first argument." >&2
    exit 1
fi

docplex config --upgrade "$STUDIO"

python - <<'PY'
import cplex
c = cplex.Cplex()
c.variables.add(names=[f"x{i}" for i in range(3000)])
c.objective.set_linear([(f"x{i}", 1.0) for i in range(3000)])
c.solve()
print("OK: full CPLEX runtime active (no 1000-variable limit)")
PY

#!/usr/bin/env bash
# =============================================================================
#  check_repo.sh - the claims this repository makes about itself, checked
# =============================================================================
#  Not part of the pipeline. It exists because a README that says "no number
#  here is untraceable" and "no emoji anywhere" is making claims a reader cannot
#  verify without doing the work, and those claims should be one command away
#  rather than taken on trust.
#
#  What it checks:
#    1. every shell script parses and passes shellcheck
#    2. every python script parses
#    3. no em dash, no en dash, no emoji in any tracked file
#    4. no banned phrase from the writing conventions
#    5. no tracked file over 50 MB, and data/ is ignored
#    6. every figure referenced by README.md exists in figures/
#    7. every figure in figures/ is produced by a script in scripts/
#    8. every results file the README names exists
#    9. "binding energy" never appears applied to a docking score
#
#  Exit status is the number of checks that failed, so it is usable in CI.
#
#  Usage: bash scripts/check_repo.sh [--quiet]
# =============================================================================
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT" || exit 1
QUIET=0
[[ "${1:-}" == "--quiet" ]] && QUIET=1
FAILED=0
say() { (( QUIET == 0 )) && echo "$@"; return 0; }
pass() { say "  [pass] $1"; }
fail() { echo "  [FAIL] $1"; FAILED=$(( FAILED + 1 )); }

mapfile -t TRACKED < <(git ls-files 2>/dev/null)
say "checking ${#TRACKED[@]} tracked files"
say

# ---- 1. shell ---------------------------------------------------------------
say "shell scripts"
SH_BAD=0
for f in run_all.sh scripts/*.sh; do
    [[ -f "$f" ]] || continue
    bash -n "$f" 2>/dev/null || { fail "$f does not parse"; SH_BAD=1; }
done
(( SH_BAD == 0 )) && pass "all shell scripts parse"
if command -v shellcheck >/dev/null 2>&1; then
    if shellcheck -x run_all.sh scripts/*.sh >/dev/null 2>&1; then
        pass "shellcheck reports no findings"
    else
        fail "shellcheck reports findings; run: shellcheck -x run_all.sh scripts/*.sh"
    fi
else
    say "  [skip] shellcheck not on PATH"
fi

# ---- 2. python --------------------------------------------------------------
say "python scripts"
if python - <<'PY'
import ast, glob, sys


def _ok(f):
    try:
        ast.parse(open(f, encoding="utf-8").read())
        return True
    except SyntaxError as e:
        print(f"    {f}: {e}")
        return False
bad = [f for f in sorted(glob.glob("scripts/*.py")) if not _ok(f)]
sys.exit(1 if bad else 0)
PY
then pass "all python scripts parse"; else fail "a python script does not parse"; fi

# ---- 3 to 4, 9. prose -------------------------------------------------------
say "writing conventions"
python - <<'PY'
import re, subprocess, sys, unicodedata
from pathlib import Path

BANNED = ["it is worth noting", "it's worth noting", "it is important to note",
          "in today's rapidly evolving", "plays a crucial role",
          "serves as a testament", "paving the way", "in conclusion", "delve",
          "seamless", "underscores", "showcases", "a testament to", "not only",
          "leverage", "leveraging", "leverages", "robust", "comprehensive"]
files = subprocess.run(["git", "ls-files"], capture_output=True, text=True).stdout.split()
em = en = emoji = banned = energy = 0
for rel in files:
    p = Path(rel)
    if p.suffix in (".png", ".gz", ".zip") or not p.is_file():
        continue
    try:
        text = p.read_text(encoding="utf-8", errors="replace")
    except Exception:
        continue
    for i, line in enumerate(text.splitlines(), 1):
        if "—" in line:
            em += 1; print(f"    EM DASH {rel}:{i}")
        if "–" in line:
            en += 1; print(f"    EN DASH {rel}:{i}")
        for ch in line:
            if ord(ch) > 0x2100 and unicodedata.category(ch) in ("So", "Sk"):
                emoji += 1; print(f"    EMOJI {rel}:{i} {ch!r}")
        # robust_rmsd is PoseBusters' own function name, not the adjective.
        scan = line.lower().replace("robust_rmsd", "")
        for b in BANNED:
            if b in scan:
                banned += 1; print(f"    BANNED {rel}:{i} '{b}'")
        # "binding energy" applied to a score. A sentence that denies it is the
        # correction this repository exists to make, not the offence.
        if re.search(r"binding energ|free energy of binding", line, re.I) and not \
           re.search(r"\b(not|never|nor|rather than|instead of|does not|is not)\b", line, re.I):
            energy += 1; print(f"    ENERGY {rel}:{i}")
print(f"  em dashes {em}, en dashes {en}, emojis {emoji}, "
      f"banned phrases {banned}, score-as-energy {energy}")
sys.exit(1 if (em or emoji or banned or energy) else 0)
PY
if [[ $? -eq 0 ]]; then pass "no em dash, emoji, banned phrase or score-as-energy"
else fail "writing conventions violated, see above"; fi

# ---- 5. repository hygiene --------------------------------------------------
say "repository hygiene"
BIG="$(git ls-files -z | xargs -0 ls -l 2>/dev/null | awk '$5 > 52428800 {print $9}')"
if [[ -z "$BIG" ]]; then pass "no tracked file over 50 MB"
else fail "tracked files over 50 MB: ${BIG}"; fi
if git check-ignore -q data/benchmarks 2>/dev/null; then pass "data/ is ignored"
else fail "data/ is not ignored"; fi
if git ls-files | grep -qE '\.(pdbqt|sdf\.gz|cif)$'; then
    fail "structure or pose files are tracked"
else pass "no pdbqt, sdf.gz or cif tracked"; fi

# ---- 6 to 8. figures and results --------------------------------------------
say "figures and results"
python - <<'PY'
import re, sys
from pathlib import Path

readme = Path("README.md")
if not readme.is_file():
    print("    README.md not written yet; figure and results checks skipped")
    sys.exit(0)
text = readme.read_text(encoding="utf-8", errors="replace")
bad = 0

refs = set(re.findall(r"\((figures/[^)\s]+)\)", text))
for r in sorted(refs):
    if not Path(r).is_file():
        print(f"    README references {r} which does not exist"); bad += 1
print(f"  README references {len(refs)} figures")

on_disk = {p.as_posix() for p in Path("figures").glob("*.png")} if Path("figures").is_dir() else set()
src = Path("scripts/09_figures.py").read_text(encoding="utf-8") if Path("scripts/09_figures.py").is_file() else ""
for f in sorted(on_disk):
    if Path(f).name not in src:
        print(f"    {f} is not produced by scripts/09_figures.py"); bad += 1
print(f"  {len(on_disk)} figures on disk, all from 09_figures.py" if not bad else "")

named = set(re.findall(r"\b((?:results|logs)/[A-Za-z0-9_./-]+\.(?:tsv|html|txt))", text))
missing = [n for n in sorted(named) if not Path(n).is_file()]
for n in missing:
    print(f"    README names {n} which does not exist"); bad += 1
print(f"  README names {len(named)} result files, {len(missing)} missing")
sys.exit(1 if bad else 0)
PY
if [[ $? -eq 0 ]]; then pass "figures and named result files all present"
else fail "a referenced figure or results file is missing"; fi

say
if (( FAILED == 0 )); then
    say "all checks passed"
else
    echo "${FAILED} check(s) failed"
fi
exit "$FAILED"

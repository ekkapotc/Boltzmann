#!/bin/sh
#
# Rename the library.  Modelled on Maxwell's nullify.sh: dry-run by default,
# bounded to this tree's own sources, and it refuses rather than guesses.
#
#     sh rename.sh boltzmann              show what would change, touch nothing
#     sh rename.sh boltzmann --apply      make the change, then build and test
#     sh rename.sh boltzmann --apply --no-verify
#
# WHAT IT TOUCHES.  inc/*.hpp, src/*.cpp, the public header,
# examples/*/main.cpp, backup_20260911/*/main.cpp, the
# Makefiles, and the .md notes.  Nothing else.  Not build/, not build-asan/.
#
# WHAT IT CHANGES, and why each one is a separate rule rather than one sed:
#
#   namespace demon            -> namespace <new>      identifier, word-bounded
#   demon::                    -> <new>::              qualified names
#   DEMON_INCLUDE, DEMON_*     -> <NEW>_*              include guards and macros
#   demon.hpp                  -> <new>.hpp            the public header, and
#                                                      every #include of it
#   demon.a / libdemon.a       -> lib<new>.a           the archive
#   -ldemon                    -> -l<new>              the linker flag
#   mpi_demon_comm             -> mpi_<new>_comm       identifier COMPONENTS:
#   mpi_started_by_demon       -> ..._by_<new>         demon glued to a
#   demon_mpi_shutdown         -> <new>_mpi_shutdown   neighbour by _ or -
#   "demon::xxx: no tape..."   -> "<new>::xxx: ..."    diagnostics the user reads
#
# THE COMPONENT RULES ARE WHY THIS IS NOT ONE EXPRESSION, and they were added
# after the first --apply proved the point.  \bdemon\b does not match inside
# -ldemon (the boundary fails after "l") or inside mpi_demon_comm (underscore
# is a word character), so the first run renamed the namespace, the header and
# the archive and left four identifiers and one linker flag behind.  The build
# failed at the link step with "cannot find -ldemon", which is the good case;
# a half-rename that still COMPILED is the one to be afraid of, and that is
# what the completeness check at the end exists to prevent.
#
# The component rules all require an adjacent _ or the -l prefix, so the word
# "demonstration" in the design note is still left alone -- verified.
#
# WHY NOT sed -i s/demon/<new>/g.  Three reasons:
#
#   1  A bare-token match.  \bdemon\b, so an identifier that merely CONTAINS
#      the letters is left alone.  This tree has none today; a substring match
#      would silently corrupt one tomorrow.
#
#   2  The word appears in PROSE as well as in code -- "Maxwell's demon is a
#      gatekeeper" in the README is about the thought experiment, not about the
#      namespace, and renaming it would turn an explanation into nonsense.
#      Lines matching PROSE_RE are reported AND excluded from every
#      substitution, by using it as a sed address.  Reporting alone is not
#      enough and was the first version's bug; see the apply block.
#
#   3  The file demon.hpp has to be renamed as well as rewritten, and git has
#      to be told, or the history of the public header is lost.
#
# Idempotent: run it twice and the second run finds nothing to do.
#
set -e

NEW="$1"
APPLY=""
VERIFY="yes"

shift 2>/dev/null || true
for a in "$@"; do
  case "$a" in
    --apply)      APPLY="yes" ;;
    --no-verify)  VERIFY="" ;;
    *) echo "rename.sh: unknown option $a" >&2; exit 2 ;;
  esac
done

if [ -z "$NEW" ]; then
  echo "usage: sh rename.sh <newname> [--apply] [--no-verify]" >&2
  exit 2
fi

case "$NEW" in
  [a-z][a-z0-9_]*) ;;
  *) echo "rename.sh: '$NEW' is not a lowercase C++ identifier" >&2; exit 2 ;;
esac

UPPER=$(echo "$NEW" | tr '[:lower:]' '[:upper:]')

# CODE AND BUILD FILES, plus README.  These are rewritten.
#
# README has no extension, so *.md does not reach it -- it was missed on the
# first dry run, which is what dry runs are for.
# backup_20260911/*/main.cpp is in here deliberately.  Those six examples do
# not compile and never will in that state, so renaming them changes nothing
# that runs -- but leaving them out means a future rename skips them AND the
# completeness check never looks at them, so the tree would quietly end up
# with two spellings of the namespace in it.  Frozen is not the same as
# invisible.
FILES=$(ls inc/*.hpp src/*.cpp demon.hpp boltzmann.hpp README Makefile \
           src/Makefile examples/Makefile examples/common.mk \
           examples/*/main.cpp examples/*/Makefile \
           backup_20260911/*/main.cpp backup_20260911/README 2>/dev/null || true)

# DESIGN NOTES.  These are reported and NEVER rewritten.
#
# A note that records the rename has to go on saying "demon" -- it is
# describing what the names used to be, and rewriting it would turn a record
# of a rename into a description of a rename that never happened.  The same
# goes for rename.sh itself, which is why it is in neither list.  Update the
# notes by hand, once, at the same time as the rename.
NOTES=$(ls *.md *.txt 2>/dev/null || true)

echo "renaming demon -> $NEW  (DEMON -> $UPPER)"
echo

# ---- rule 2: prose that must not be rewritten ----------------------------
#
# One pattern, used in three places: to report the prose, to keep it out of
# the "what would change" count, and to keep it out of the completeness check
# at the end.  When the README gained a paragraph explaining the name -- which
# says "Maxwell's demon is one gatekeeper" and must keep saying it -- the
# completeness check flagged the tree as half-renamed.  Three copies of the
# same idea in three slightly different regexes is how that happens, so there
# is one copy.
#
PROSE_RE="maxwell'?s demon|demon is a gatekeeper|the demon\b|demonstrat"

PROSE=$(grep -n -i -E "$PROSE_RE" $FILES 2>/dev/null || true)

if [ -n "$PROSE" ]; then
  echo "PROSE in code/build files, left alone -- check these by hand:"
  echo "$PROSE" | sed 's/^/    /'
  echo
fi

if [ -n "$NOTES" ]; then
  NOTEHITS=$(grep -H -c -i "demon" $NOTES 2>/dev/null | grep -v ':0$' || true)
  if [ -n "$NOTEHITS" ]; then
    echo "DESIGN NOTES, never rewritten -- update these by hand:"
    echo "$NOTEHITS" | sed 's/^/    /'
    echo
  fi
fi

# ---- what would change ---------------------------------------------------
#
# Detect with the SAME breadth the rules cover, not with -w.  Using -w here
# was the second half of the -ldemon bug: after the word-bounded rules had
# run, -w found nothing left, the script said "nothing to do" and exited 0 --
# reporting success on a tree that no longer linked.  The detector and the
# rules have to agree about what counts as an occurrence.
#
HITS=$(grep -H -c -i "demon" $FILES 2>/dev/null | grep -v ':0$' || true)

# ...minus the prose, which is about the thought experiment and stays.
# The `&& echo` below is the last command in the loop body, so on a file whose
# only matches are prose it returns false -- and under `set -e` a command
# substitution ending in a false command kills the script.  It did: the tree
# was fully renamed, the check found nothing left, and rename.sh exited 1
# without printing anything.  Hence the explicit `|| true` per iteration.
HITS=$(for h in $HITS; do
         f=${h%%:*}
         n=$(grep -i "demon" "$f" 2>/dev/null | grep -v -c -i -E "$PROSE_RE" || true)
         if [ "${n:-0}" -gt 0 ]; then echo "$f:$n"; fi
       done || true)

if [ -z "$HITS" ]; then
  echo "nothing to do."
  exit 0
fi

echo "files with occurrences:"
echo "$HITS" | sed 's/^/    /'
echo

if [ -z "$APPLY" ]; then
  echo "dry run.  re-run with --apply to make the change."
  exit 0
fi

# ---- apply ---------------------------------------------------------------
#
# EVERY SUBSTITUTION IS GUARDED BY THE PROSE ADDRESS.  It has to be, and the
# first version was not: it printed "PROSE, left alone" and then rewrote the
# prose anyway, because reporting a line and excluding it from the sed are two
# different things and only the first was implemented.  The bug was invisible
# in practice -- the README paragraph explaining the name was written after
# the last --apply, so the sed never ran over it -- and it surfaced only when
# the pre-rename tree was reconstructed and the rename re-run, which turned
# "Maxwell's demon" into "Maxwell's boltzmann".  A round trip is worth more
# than an inspection.
#
# /re/I! is "on lines NOT matching re, case-insensitively", so a prose line is
# passed through untouched no matter which rule would otherwise have hit it.
#
for f in $FILES; do
  [ -f "$f" ] || continue
  sed -E -i "/${PROSE_RE}/I!{
    s/\bDEMON_/${UPPER}_/g
    s/\bdemon\.hpp\b/${NEW}.hpp/g
    s/\blibdemon\b/lib${NEW}/g
    s/\bdemon\.a\b/${NEW}.a/g
    s/-ldemon\b/-l${NEW}/g
    s/_demon_/_${NEW}_/g
    s/_demon\b/_${NEW}/g
    s/\bdemon_/${NEW}_/g
    s/\bdemon\b/${NEW}/g
  }" "$f"
done

if [ -f demon.hpp ]; then
  if [ -d .git ]; then git mv demon.hpp "${NEW}.hpp"; else mv demon.hpp "${NEW}.hpp"; fi
fi

echo "renamed.  public header is now ${NEW}.hpp, archive lib${NEW}.a"

# ---- completeness check --------------------------------------------------
#
# A half-rename that still links is the dangerous outcome, so look for
# anything left over rather than trusting the rules.  Prose is excluded by the
# same pattern used above; anything else is reported and the script exits
# non-zero.
#
LEFT=$(grep -rn -i "demon" $FILES 2>/dev/null | grep -v -i -E "$PROSE_RE" || true)

if [ -n "$LEFT" ]; then
  echo
  echo "INCOMPLETE -- these still mention demon:"
  echo "$LEFT" | sed 's/^/    /'
  echo
  echo "add a rule for each form above and re-run; the script is idempotent."
  exit 1
fi

echo "completeness check: no occurrence of 'demon' remains."

if [ -n "$VERIFY" ]; then
  echo
  echo "verifying: make clean && make && make ranks"
  make clean >/dev/null
  make
  make ranks
fi

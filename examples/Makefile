#
# Every subdirectory that contains a Makefile and no SKIP file is an example.
# Drop a new directory in with a main.cpp and a Makefile that includes
# ../common.mk, and `make`, `make clean` and `make test` pick it up with no
# edit here.
#
# This replaces a hand-maintained list that had drifted: it named eleven
# directories where the tree held seventeen.  The six it left out no longer
# compile -- they call an API removed years ago -- and they are frozen in
# ../backup_20260911/ rather than skipped in place, which is what Maxwell
# does with its own pre-checkpoint examples.
#
# SKIP is kept as a facility even though nothing currently uses it: a
# directory with a SKIP file in it is passed over, and the first line of that
# file is printed by `make list`.  It is the right answer for an example that
# is being worked on rather than one that has been retired, and it costs two
# lines.
#
ALL_DIRS  := $(sort $(patsubst %/,%,$(dir $(wildcard */Makefile))))
SKIPPED   := $(sort $(patsubst %/,%,$(dir $(wildcard */SKIP))))
EXAMPLES_ALL := $(filter-out $(SKIPPED),$(ALL_DIRS))

# The assertion suites are seconds where the rest are minutes, so they go
# first.  invariants leads: it is the one that says whether the library still
# computes the same thing, so a failure there explains the rest.
FAST      := invariants fdcheck regress
EXAMPLES  := $(foreach f,$(FAST),$(filter $(f),$(EXAMPLES_ALL))) $(filter-out $(FAST),$(EXAMPLES_ALL))

all clean:
	@for d in $(EXAMPLES); do $(MAKE) --no-print-directory -C $$d $@ || exit 1; done

test: all
	@for d in $(EXAMPLES); do \
	  echo "=================================================== $$d"; \
	  $(MAKE) --no-print-directory -C $$d run || exit 1; \
	done
	@echo "=================================================== all examples passed"

list:
	@for d in $(EXAMPLES); do \
	  printf "  %-12s %s\n" "$$d" "`$(MAKE) --no-print-directory -C $$d echo-args`"; \
	done
	@for d in $(SKIPPED); do printf "  %-12s SKIPPED: %s\n" "$$d" "`head -1 $$d/SKIP`"; done

.PHONY: all clean test list

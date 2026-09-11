#
# A pipelined-parallel automatic differentiation library: Jacobians by sparse
# graph vertex elimination, with the tape partitioned across MPI ranks.
#
#   inc/         headers.  The public surface is boltzmann.hpp, which pulls in
#                Typedefs, Active, Jacobian, API and BreakException.  Vertex,
#                Edge, Process and the two arenas are internals and are
#                deliberately not reachable from it -- nor is <mpi.h>.
#   src/         the library.  Six translation units, built by src/Makefile
#                into build/lib/libboltzmann.a.
#   examples/    one directory per example: a main.cpp and a Makefile of one
#                to four lines that includes ../common.mk.  examples/Makefile
#                discovers them from the filesystem -- there is no list to
#                edit when one is added or removed.  A directory containing
#                a file called SKIP is passed over.
#   build/       everything generated.  Nothing is written into the source
#                tree, and `make clean` is `rm -rf build`.
#   backup_20260911/
#                six examples that predate the checkpointing protocol and no
#                longer compile, frozen with their own README.  Not built,
#                not discovered, kept for the record -- as Maxwell keeps its
#                own pre-checkpoint examples.
#
# Targets
#   make            the library and every example
#   make test       build, then run the whole suite at NP ranks
#   make list       every example discovered, with the arguments it runs with
#   make ranks      run the assertion suites at 1, 2, 3 and 4 ranks.  THE
#                   IMPORTANT ONE: the Jacobian must not depend on how many
#                   ranks computed it, and nothing else in the suite checks
#                   that, because every other target runs at one NP.
#   make sanitize   rebuild under AddressSanitizer in a SEPARATE tree and run
#                   the assertion suites
#   make clean      remove build/
#
# Ported from Maxwell's build tree, with one addition: NP.  Every example is
# launched through $(MPIRUN) -np $(NP), so the whole suite can be moved to a
# different rank count from the command line -- make test NP=8.
#
BUILD_DIR := $(CURDIR)/build
export BUILD_DIR

# make PREDEFINES CXX as "g++", so "CXX ?= mpic++" silently does nothing --
# ?= only fires when the variable is undefined, and a built-in default counts
# as defined.  This is the same trap Maxwell hit with CC.  Override the
# built-in default only, so `make CXX=mpiicpc` and an exported CXX both win.
ifeq ($(origin CXX),default)
CXX := mpic++
endif
MPIRUN    ?= mpirun
NP        ?= 2
export CXX MPIRUN NP

STD       ?= -std=gnu++0x
WARN      ?= -Wall -Wextra
LIB_FLAGS ?= -g $(STD) $(WARN)

.PHONY: all lib examples test list ranks sanitize clean help

all: examples

lib:
	@$(MAKE) --no-print-directory -C src FLAGS="$(LIB_FLAGS)"

examples: lib
	@$(MAKE) --no-print-directory -C examples all

test: all
	@$(MAKE) --no-print-directory -C examples test

list:
	@$(MAKE) --no-print-directory -C examples list

# --------------------------------------------------------------- ranks ----
#
# The invariant that only a parallel library has: the Jacobian a program
# computes must not depend on the number of ranks that computed it, or on the
# memory budget, or on the two together.  `make test` runs at one NP and
# cannot see a violation; this runs the assertion suites across four.
#
RANK_EXAMPLES ?= invariants fdcheck regress
RANK_COUNTS   ?= 1 2 3 4

ranks: all
	@for np in $(RANK_COUNTS); do \
	  for d in $(RANK_EXAMPLES); do \
	    echo "=================================================== $$d  (np=$$np)"; \
	    $(MAKE) --no-print-directory -C examples/$$d run NP=$$np || exit 1; \
	  done; \
	done
	@echo "=================================================== identical at $(RANK_COUNTS) ranks"

# ------------------------------------------------------------ sanitize ----
#
# A separate BUILD_DIR rather than a flags stamp: make compares timestamps
# only, so sanitizing into build/ leaves every binary newer than its source
# and a following plain `make` rebuilds nothing -- you keep running
# instrumented code without being told.
#
SAN_FLAGS    ?= -fsanitize=address -fno-omit-frame-pointer
SAN_BUILD    ?= $(CURDIR)/build-asan
SAN_EXAMPLES ?= invariants regress fdcheck

sanitize:
	@$(MAKE) --no-print-directory -C src \
	    BUILD_DIR="$(SAN_BUILD)" \
	    FLAGS="-g $(STD) $(WARN) -O1 $(SAN_FLAGS)"
	@for d in $(SAN_EXAMPLES); do \
	  $(MAKE) --no-print-directory -C examples/$$d all \
	      BUILD_DIR="$(SAN_BUILD)" \
	      BASE_FLAGS="-g $(STD) $(WARN) -O1 $(SAN_FLAGS)" || exit 1; \
	done
	@for d in $(SAN_EXAMPLES); do \
	  echo "=================================================== $$d (asan)"; \
	  $(MAKE) --no-print-directory -C examples/$$d run BUILD_DIR="$(SAN_BUILD)" || exit 1; \
	done
	@echo "=================================================== sanitizer clean"

clean:
	@$(MAKE) --no-print-directory -C examples clean
	@rm -rf $(BUILD_DIR) $(SAN_BUILD)

help:
	@echo "make            build build/lib/libboltzmann.a and every example into build/bin"
	@echo "make test       build, then run the whole suite at NP=$(NP) ranks"
	@echo "make ranks      run the assertion suites at $(RANK_COUNTS) ranks"
	@echo "make list       list the examples and the arguments they run with"
	@echo "make sanitize   rebuild into build-asan/ under ASan and run: $(SAN_EXAMPLES)"
	@echo "make clean      remove build/ and build-asan/"

#
# Shared build rules for every example.  An example directory needs nothing
# but main.cpp and a Makefile of one to four lines:
#
#     RUNARGS     = 20 1          # what `make run` passes to the binary
#     NP          = 4             # rank count, if this example needs a fixed one
#     EXTRA_FLAGS = -O2           # optional
#     include ../common.mk
#
# ../Makefile discovers the directory automatically.  Drop a file called SKIP
# in a directory to have the suite pass over it.
#
# Nothing is written into the example directory.  The object goes to
# $(BUILD_DIR)/obj/examples/<name>.o and the binary to $(BUILD_DIR)/bin/<name>,
# named after the directory, so the binaries can share one bin/ instead of
# all being called "main".
#
NAME       := $(notdir $(CURDIR))
ROOT_DIR   := $(abspath ../..)

BUILD_DIR  ?= $(ROOT_DIR)/build
OBJ_DIR     = $(BUILD_DIR)/obj/examples
BIN_DIR     = $(BUILD_DIR)/bin
LIB_DIR     = $(BUILD_DIR)/lib
LIBRARY     = $(LIB_DIR)/libboltzmann.a

# An example's object depended on its main.cpp and on the archive, so editing
# boltzmann.hpp or a header under inc/ rebuilt the example only if the library
# happened to be re-archived as well.
API_HEADERS = $(ROOT_DIR)/boltzmann.hpp $(wildcard $(ROOT_DIR)/inc/*.hpp)

PROG        = $(BIN_DIR)/$(NAME)
OBJ         = $(OBJ_DIR)/$(NAME).o

# CXX, not CC: make predefines CC as "cc", so "CC ?= mpic++" silently does
# nothing and the examples link with the C driver -- no libstdc++, and a page
# of undefined references.
# make PREDEFINES CXX as "g++", so "CXX ?= mpic++" silently does nothing --
# ?= only fires when the variable is undefined, and a built-in default counts
# as defined.  This is the same trap Maxwell hit with CC.  Override the
# built-in default only, so `make CXX=mpiicpc` and an exported CXX both win.
ifeq ($(origin CXX),default)
CXX := mpic++
endif
MPIRUN     ?= mpirun
NP         ?= 2

# Set and exported by ../../Makefile, which explains what it is for.  The
# fallback is here so that `make run` works in an example directory on its own,
# without going through the top-level Makefile: Open MPI needs --oversubscribe
# to accept -np greater than the core count, MPICH rejects the flag.
MPIRUN_FLAGS ?= $(shell $(MPIRUN) --version 2>&1 | grep -qi 'open[ -]*mpi' && echo --oversubscribe)
BASE_FLAGS ?= -g -std=gnu++0x -Wall -Wextra
EXTRA_FLAGS ?=
FLAGS       = $(BASE_FLAGS) $(EXTRA_FLAGS)
LIBS        = -L$(LIB_DIR) -lboltzmann
RUNARGS    ?=

all: $(PROG)

$(PROG): $(OBJ) $(LIBRARY) | $(BIN_DIR)
	@echo "  LD  $@"
	@$(CXX) $(FLAGS) -o $@ $(OBJ) $(LIBS)

$(OBJ): main.cpp $(API_HEADERS) $(LIBRARY) | $(OBJ_DIR)
	@echo "  CXX $(NAME)/main.cpp"
	@$(CXX) $(FLAGS) -c $< -o $@

# So that a single example builds from a cold tree without going via the top.
$(LIBRARY):
	@$(MAKE) --no-print-directory -C $(ROOT_DIR)/src

$(OBJ_DIR) $(BIN_DIR):
	@mkdir -p $@

run: all
	$(MPIRUN) $(MPIRUN_FLAGS) -np $(NP) $(PROG) $(RUNARGS)

echo-args:
	@echo "$(MPIRUN_FLAGS) -np $(NP) $(RUNARGS)"

clean:
	@rm -rf $(OBJ) $(PROG) test.*.er output*.txt job.sh core* *.btr *.o main

.PHONY: all clean run echo-args

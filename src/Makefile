#
# The library.  Objects and the archive go under $(BUILD_DIR); nothing is
# written into src/.
#
# BUILD_DIR is exported by the top-level Makefile; the fallback keeps a bare
# `cd src && make` working.
#
ROOT_DIR  := $(abspath ..)
BUILD_DIR ?= $(ROOT_DIR)/build
OBJ_DIR    = $(BUILD_DIR)/obj/src
LIB_DIR    = $(BUILD_DIR)/lib
LIBRARY    = $(LIB_DIR)/libboltzmann.a

# make PREDEFINES CXX as "g++", so "CXX ?= mpic++" silently does nothing --
# ?= only fires when the variable is undefined, and a built-in default counts
# as defined.  This is the same trap Maxwell hit with CC.  Override the
# built-in default only, so `make CXX=mpiicpc` and an exported CXX both win.
ifeq ($(origin CXX),default)
CXX := mpic++
endif
CC        := $(CXX)
FLAGS     ?= -g -std=gnu++0x -Wall -Wextra
INC_DIR    = $(ROOT_DIR)/inc

# Every object depends on EVERY header.  The old rule named one or two
# headers each, so editing Typedefs.hpp -- which all six include -- rebuilt
# almost nothing, and editing a header that a given .cpp did not name left
# that object stale against a changed struct layout.  That is a silent ODR
# violation, and in this library it means one translation unit thinking a
# Vertex is 120 bytes while another thinks it is 72.  Six files do not need
# fine-grained dependencies.
HEADERS    = $(wildcard $(INC_DIR)/*.hpp) $(ROOT_DIR)/boltzmann.hpp

NAMES      = API Active Vertex Edge Process BreakException
OBJECTS    = $(addprefix $(OBJ_DIR)/,$(addsuffix .o,$(NAMES)))

all: $(LIBRARY)

$(LIBRARY): $(OBJECTS) | $(LIB_DIR)
	@echo "  AR  $@"
	@rm -f $@
	@ar -cr $@ $(OBJECTS)

$(OBJ_DIR)/%.o: %.cpp $(HEADERS) | $(OBJ_DIR)
	@echo "  CXX src/$<"
	@$(CC) $(FLAGS) -c $< -I$(INC_DIR) -o $@

$(OBJ_DIR) $(LIB_DIR):
	@mkdir -p $@

clean:
	@rm -rf $(OBJ_DIR) $(LIBRARY) *.o *.a

.PHONY: all clean

################################################################################

# Sources and targets
# use TARGET for host only programs (no GPU)
# use NVTARGET for GPU programs
TARGET = solution
# NVTARGET = solution
MODULES = $(if $(wildcard solution.*),solution,template)
OBJECTS = $(addsuffix .o,$(MODULES))

################################################################################

include ../common.mak
OPT += -fopenmp
LDFLAGS += -fopenmp

# Include "data" as an order-only prerequisite to generate data
# e.g. run: all | data
.PHONY: run
run: all
	for t in 1 2 {4..20..4}; do \
	  { echo "Max threads: $$t"; \
	  OMP_NUM_THREADS=$$t ./solution;} \
	done

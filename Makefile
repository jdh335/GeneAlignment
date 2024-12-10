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

.PHONY: cuda
cuda: all
	$(MAKE) -C ./CUDA run

.PHONY: run
run: all
	rm -f run_log_omp.txt
	for t in 1 2 {4..20..4}; do \
	  { echo "Max threads: $$t"; \
	  OMP_NUM_THREADS=$$t ./solution;} \
		| tee -a run_log_omp.txt; \
	done
	make cuda
	

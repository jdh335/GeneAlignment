################################################################################

# Sources and targets
# use TARGET for host only programs (no GPU)
# use NVTARGET for GPU programs
# TARGET = solution
NVTARGET = solution
MODULES = $(if $(wildcard solution.*),solution,template)
OBJECTS = $(addsuffix .o,$(MODULES))

################################################################################

include ../../common.mak
HTK_DIR = ../../libhtk

.PHONY: run
run: all 
	  ./solution \
	done

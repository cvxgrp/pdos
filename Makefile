include pdos.mk
# If you wish to use this experimental build, you will have to provide
# the path to Pardiso and a BLAS library.
PATH_TO_PARDISO=/home/echu/pardiso
PATH_TO_BLAS=/home/echu/OpenBLAS
LDFLAGS += $(PATH_TO_PARDISO)/libpardiso412-GNU450-X86-64.so -fopenmp
# to link to system blas and lapack uncomment the following line
# LDFLAGS += -llapack -lblas

# to link to OpenBLAS uncomment the following line 
LDFLAGS += $(PATH_TO_BLAS)/libopenblas.so

OBJECTS = pdos.o util.o cones.o cs.o
AMD_SOURCE = $(wildcard direct/amd_*.c)
DIRECT_OBJECTS = direct/ldl.o $(AMD_SOURCE:.c=.o)
TARGETS = demo_direct demo_indirect

.PHONY: default
default: lib/libpdos.so bin/pdos

pdos.o 		: pdos.h globals.h common.h linAlg.h
util.o		: util.h pdos.h
cones.o		: cones.h globals.h
cs.o		: cs.h globals.h

direct/private.o	: direct/private.h common.h
direct/ldl.o		: direct/ldl.h
direct/amd_1.o		: direct/amd_internal.h direct/amd.h
direct/amd_2.o		: direct/amd_internal.h direct/amd.h
direct/amd_aat.o	: direct/amd_internal.h direct/amd.h
direct/amd_control.o	: direct/amd_internal.h direct/amd.h
direct/amd_defaults.o 	: direct/amd_internal.h direct/amd.h
direct/amd_dump.o	: direct/amd_internal.h direct/amd.h
direct/amd_global.o	: direct/amd_internal.h direct/amd.h
direct/amd_info.o	: direct/amd_internal.h direct/amd.h
direct/amd_order.o	: direct/amd_internal.h direct/amd.h
direct/amd_post_tree.o	: direct/amd_internal.h direct/amd.h
direct/amd_postorder.o	: direct/amd_internal.h direct/amd.h
direct/amd_preprocess.o	: direct/amd_internal.h direct/amd.h
direct/amd_valid.o	: direct/amd_internal.h direct/amd.h

indirect/private.o	: indirect/private.h common.h

lib/libpdos.so: $(OBJECTS) direct/private.o  $(DIRECT_OBJECTS)
	mkdir -p lib
	$(ARCHIVE) lib/libpdosdir.a $^
	- $(RANLIB) lib/libpdosdir.a

bin/pdos: run_pdos.c lib/libpdosdir.a
	mkdir -p bin
	$(CC) $(CFLAGS) -DDEMO_PATH="\"$(CURDIR)/data/portfolio_test_01_10x100\"" -o $@ $^ $(LDFLAGS)

.PHONY: clean purge
clean:
	@rm -rf $(TARGETS) $(OBJECTS) $(DIRECT_OBJECTS) direct/private.o indirect/private.o

purge: clean
	@rm -rf bin lib

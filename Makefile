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
TARGETS = pdos 

.PHONY: default
default: lib/libpdos.a bin/pdos

pdos.o 		: pdos.h globals.h common.h linAlg.h
util.o		: util.h pdos.h
cones.o		: cones.h globals.h
cs.o		  : cs.h globals.h

pardiso/private.o	    : pardiso/private.h common.h

lib/libpdos.a: $(OBJECTS) pardiso/private.o  $(DIRECT_OBJECTS)
	mkdir -p lib
	$(ARCHIVE) lib/libpdos.a $^
	- $(RANLIB) lib/libpdos.a

bin/pdos: run_pdos.c lib/libpdos.a
	mkdir -p bin
	$(CC) $(CFLAGS) -DDEMO_PATH="\"$(CURDIR)/data/portfolio_test_01_10x100\"" -o $@ $^ $(LDFLAGS)

.PHONY: clean purge
clean:
	@rm -rf $(TARGETS) $(OBJECTS) pardiso/private.o

purge: clean
	@rm -rf bin lib

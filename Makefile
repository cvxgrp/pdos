include pdos.mk

# libraries for AMD and LDL, needed for direct method:
EINCS = -Idirect/external/LDL/Include -Idirect/external/AMD/Include -Idirect/external/SuiteSparse_config 
AMD_LIB = direct/external/AMD/Lib/libamd.a 
LDL_LIB = direct/external/LDL/Lib/libldl.a

OBJECTS = pdos.o linAlg.o common.o util.o cones.o cs.o
TARGETS = demo_direct demo_indirect

default: amd ldl pdos_direct pdos_indirect demo_direct demo_indirect

packages: amd ldl

amd:
	(cd direct/external/AMD ; $(MAKE))

ldl:
	(cd direct/external/LDL ; $(MAKE))

pdos.o 		: pdos.h
common.o	: common.h
linAlg.o	: linAlg.h
util.o		: util.h
cones.o		: cones.h
cs.o			: cs.h

direct/private.o: direct/private.c direct/private.h
	$(CC) $(CFLAGS) $(EINCS) -c -o $@ direct/private.c

indirect/private.o: indirect/private.c indirect/private.h
	$(CC) $(CFLAGS) -c -o $@ indirect/private.c

pdos_direct: $(OBJECTS) direct/private.o $(AMD_LIB) $(LDL_LIB)
	mkdir -p lib
	$(ARCHIVE) lib/libpdosdir.a $^
	- $(RANLIB) lib/libpdosdir.a

pdos_indirect: $(OBJECTS) indirect/private.o
	mkdir -p lib
	$(ARCHIVE) lib/libpdosindir.a $^
	- $(RANLIB) lib/libpdosindir.a

demo_direct: run_pdos.c 
	mkdir -p bin
	$(CC) $(CFLAGS) -DDEMO_PATH="\"$(CURDIR)/data_pdos\"" -o bin/$@ $^ lib/libpdosdir.a $(LDFLAGS) 

demo_indirect: run_pdos.c 
	mkdir -p bin
	$(CC) $(CFLAGS) -DDEMO_PATH="\"$(CURDIR)/data_pdos\"" -o bin/$@ $^ lib/libpdosindir.a $(LDFLAGS) 

.PHONY: clean purge

clean:
	( cd direct/external/LDL      ; $(MAKE) clean )
	( cd direct/external/AMD      ; $(MAKE) clean )
	@rm -rf $(TARGETS) $(OBJECTS) direct/private.o indirect/private.o core Makefile.dependencies *.o *.a

purge: clean
	( cd direct/external/LDL    ; $(MAKE) purge )
	( cd direct/external/AMD    ; $(MAKE) purge )   
	@rm -rf bin lib

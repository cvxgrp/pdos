include pdos.mk

OBJECTS = pdos.o util.o cones.o cs.o
AMD_SOURCE = $(wildcard direct/amd_*.c)
DIRECT_OBJECTS = direct/ldl.o $(AMD_SOURCE:.c=.o) 
TARGETS = demo_direct demo_indirect

.PHONY: default 
	#lib/libpdosdir.a lib/libpdosindir.a bin/demo_direct bin/demo_indirect
default: lib/libpdosdir.a lib/libpdosindir.a bin/demo_direct bin/demo_indirect

pdos.o 		: pdos.h numericDef.h common.h linAlg.h
util.o		: util.h pdos.h
cones.o		: cones.h numericDef.h
cs.o			: cs.h numericDef.h

direct/private.o				: direct/private.h common.h
direct/ldl.o						: direct/ldl.h
direct/amd_1.o					: direct/amd_internal.h direct/amd.h
direct/amd_2.o					: direct/amd_internal.h direct/amd.h
direct/amd_aat.o				: direct/amd_internal.h direct/amd.h
direct/amd_control.o		: direct/amd_internal.h direct/amd.h
direct/amd_defaults.o 	: direct/amd_internal.h direct/amd.h
direct/amd_dump.o				: direct/amd_internal.h direct/amd.h
direct/amd_global.o			: direct/amd_internal.h direct/amd.h
direct/amd_info.o				: direct/amd_internal.h direct/amd.h
direct/amd_order.o			: direct/amd_internal.h direct/amd.h
direct/amd_post_tree.o	: direct/amd_internal.h direct/amd.h
direct/amd_postorder.o	: direct/amd_internal.h direct/amd.h
direct/amd_preprocess.o	: direct/amd_internal.h direct/amd.h
direct/amd_valid.o			: direct/amd_internal.h direct/amd.h

indirect/private.o	: indirect/private.h common.h

lib/libpdosdir.a: $(OBJECTS) direct/private.o  $(DIRECT_OBJECTS)
	mkdir -p lib
	$(ARCHIVE) lib/libpdosdir.a $^
	- $(RANLIB) lib/libpdosdir.a

lib/libpdosindir.a: $(OBJECTS) indirect/private.o
	mkdir -p lib
	$(ARCHIVE) lib/libpdosindir.a $^
	- $(RANLIB) lib/libpdosindir.a

bin/demo_direct: run_pdos.c lib/libpdosdir.a
	mkdir -p bin
	$(CC) $(CFLAGS) -DDEMO_PATH="\"$(CURDIR)/data_pdos\"" -o $@ $^ $(LDFLAGS) 

bin/demo_indirect: run_pdos.c lib/libpdosindir.a
	mkdir -p bin
	$(CC) $(CFLAGS) -DDEMO_PATH="\"$(CURDIR)/data_pdos\"" -o $@ $^ $(LDFLAGS) 

.PHONY: clean purge
clean:
	@rm -rf $(TARGETS) $(OBJECTS) $(DIRECT_OBJECTS) direct/private.o indirect/private.o

purge: clean 
	@rm -rf bin lib

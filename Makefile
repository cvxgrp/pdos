include pdos.mk

OBJECTS = pdos.o linAlg.o common.o util.o cones.o cs.o
AMD_SOURCE = $(wildcard direct/amd_*.c)
DIRECT_OBJECTS = direct/ldl.o $(AMD_SOURCE:.c=.o) 
TARGETS = demo_direct demo_indirect

default: pdos_direct pdos_indirect demo_direct demo_indirect

%.o : %.h

pdos_direct: $(OBJECTS) direct/private.o  $(DIRECT_OBJECTS)
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
	@rm -rf $(TARGETS) $(OBJECTS) $(DIRECT_OBJECTS) direct/private.o indirect/private.o

purge: clean 
	@rm -rf bin lib

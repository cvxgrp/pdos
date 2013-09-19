UNAME = $(shell uname -s)
CC = gcc	# use GCC for openmp support
CFLAGS = -pg -Wall -pedantic -O3 -I. -DDLONG
LDFLAGS = -lm 

ifeq ($(UNAME), Darwin)
	CFLAGS   += -std=c99
else
	CFLAGS   += -std=gnu99
endif

ifeq ($(UNAME), Linux)
# we're on a linux system, use accurate timer provided by clock_gettime()
LDFLAGS += -lrt
endif


AR = ar
ARFLAGS = rv
ARCHIVE = $(AR) $(ARFLAGS)
RANLIB = ranlib

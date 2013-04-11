UNAME = $(shell uname -s)
CC = cc
CFLAGS = -g -Wall -pedantic -O3 -I. -DDLONG -DLDL_LONG
LDFLAGS = -lm

ifeq ($(UNAME), Darwin)
	CFLAGS   += -std=c99
else
	CFLAGS   += -std=gnu99
endif

AR = ar
ARFLAGS = rv
ARCHIVE = $(AR) $(ARFLAGS)
RANLIB = ranlib

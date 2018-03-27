CC = gcc
OPTS = -Wall -Wextra -O2 -g -std=c99

ifdef HTSLIB
L=$(HTSLIB)/libhts.a -L$(HTSLIB)
I=-I$(HTSLIB)
endif


.PHONY:	all clean

.SUFFIXES: .c .o

.c.o: Makefile
	$(CC) -c $(OPTS) $< -o $@ $I

all: strand_cov

strand_cov: strand_cov.o
	$(CC) $(OPTS) $(INCLUDES) -o strand_cov strand_cov.o $L -lhts -lz -lpthread -lm -llzma -lbz2

clean:
	rm -f *.o strand_cov

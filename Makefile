INCLUDES = -I../htslib
CC = gcc
OPTS = -Wall -Wextra -O2 -g -std=c99

.PHONY:	all htslib clean

.SUFFIXES: .c .o

.c.o: Makefile
	$(CC) -c $(OPTS) $(INCLUDES) $< -o $@

all: strand_cov

htslib:
	$(MAKE) -j -C ../htslib

strand_cov: htslib strand_cov.o
	$(CC) $(OPTS) $(INCLUDES) -o strand_cov strand_cov.o ../htslib/libhts.a -lz -lpthread -lm -llzma -lbz2

clean:
	rm -f *.o strand_cov

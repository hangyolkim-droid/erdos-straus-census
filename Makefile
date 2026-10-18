CC      ?= cc
CFLAGS  ?= -O3 -march=native -Wall
LDLIBS  ?= -lm

all: c27 es_big bench

c27: src/c27.c
	$(CC) $(CFLAGS) -pthread -o $@ $< $(LDLIBS)

es_big: src/es_big.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

bench: src/bench.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f c27 es_big bench

.PHONY: all clean

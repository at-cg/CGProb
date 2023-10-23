CC=g++
CFLAGS=-O3
STD=-std=c++2b
LIBS=-lgmp -lgmpxx
OMP=-fopenmp

.PHONY: all clean

all: compute

compute: compute.cpp
	$(CC) $(CFLAGS) $(OMP) $(STD) -o compute compute.cpp $(LIBS)

clean:
	rm -f compute
import sys

# Check if correct number of arguments were passed
if len(sys.argv)!=3:
    print("Usage: python script.py path_to_include path_to_lib")
    sys.exit(1)

path_to_include = sys.argv[1]
path_to_lib = sys.argv[2]

print("CC=g++")
print("CFLAGS=-O3")
print("STD=-std=c++2b")
print("LIBS=-lgmp -lgmpxx")
print("OMP=-fopenmp")

print(".PHONY: all clean\n")

print("all: compute\n")

print("compute: compute.cpp")
print(f"\t$(CC) $(CFLAGS) $(OMP) $(STD) -o compute compute.cpp $(LIBS) -l {path_to_include} -L {path_to_lib}\n")

print("clean:")
print("\trm -f compute")

# Amaunet: High-order Lattice Perturbation Theory
#          for Non-relativistic Quantum Matter
#
# High-order Perturbation Theory Analytics
# Weak-coupling Expansion for Fermionic Contact Interactions
#
# Primary Makefile
#   - Main executable (target: all)
#   - Unit testing (target: ut)
#
# Andrew C. Loheac, Joaquin E. Drut
# Department of Physics and Astronomy
# University of North Carolina at Chapel Hill

CC=g++
CFLAGS=-g

all: amaunet

amaunet: main.o PTSymbolicObjects.o
	$(CC) $(CFLAGS) main.o PTSymbolicObjects.o -o amaunet
	
main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp
	
PTSymbolicObjects.o: PTSymbolicObjects.cpp
	$(CC) $(CFLAGS) -c PTSymbolicObjects.cpp
	
ut: unittst

unittst: UnitTesting.o PTSymbolicObjects.o
	$(CC) $(CFLAGS) UnitTesting.o PTSymbolicObjects.o -o unittst
	
UnitTesting.o: UnitTesting.cpp
	$(CC) $(CFLAGS) -c UnitTesting.cpp

clean:
	rm *o amaunet
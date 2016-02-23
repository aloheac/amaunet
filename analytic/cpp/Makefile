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
CFLAGS=-g -std=c++11

all: amaunet

amaunet: main.o PTSymbolicObjects.o PathIntegration.o
	$(CC) $(CFLAGS) main.o PTSymbolicObjects.o PathIntegration.o -o amaunet
	
main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp
	
PTSymbolicObjects.o: PTSymbolicObjects.cpp
	$(CC) $(CFLAGS) -c PTSymbolicObjects.cpp
	
PathIntegration.o: PathIntegration.cpp
	$(CC) $(CFLAGS) -c PathIntegration.cpp
	
ut: unittst

unittst: UnitTesting.o PTSymbolicObjects.o PathIntegration.o
	$(CC) $(CFLAGS) UnitTesting.o PTSymbolicObjects.o PathIntegration.o -o unittst
	
UnitTesting.o: UnitTesting.cpp
	$(CC) $(CFLAGS) -c UnitTesting.cpp

clean:
	rm *o amaunet unittst
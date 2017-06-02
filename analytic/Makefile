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
CFLAGS=-g -std=c++11 -fopenmp
LIBBOOST=-lboost_serialization

all: amaunet

amaunet: main.o PTSymbolicObjects.o PathIntegration.o Multithreading.o Debugging.o FeynmanDiagram.o ExpressionSerialization.o
	$(CC) $(CFLAGS) main.o PTSymbolicObjects.o PathIntegration.o Multithreading.o Debugging.o FeynmanDiagram.o ExpressionSerialization.o -o amaunet $(LIBBOOST)
	
main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp
	
PTSymbolicObjects.o: PTSymbolicObjects.cpp
	$(CC) $(CFLAGS) -c PTSymbolicObjects.cpp
	
PathIntegration.o: PathIntegration.cpp
	$(CC) $(CFLAGS) -c PathIntegration.cpp
	
Multithreading.o: Multithreading.cpp
	$(CC) $(CFLAGS) -c Multithreading.cpp
	
Debugging.o: Debugging.cpp
	$(CC) $(CFLAGS) -c Debugging.cpp

FeynmanDiagram.o: FeynmanDiagram.cpp
	$(CC) $(CFLAGS) -c FeynmanDiagram.cpp

ExpressionSerialization.o: ExpressionSerialization.cpp
	$(CC) $(CFLAGS) -c ExpressionSerialization.cpp
	
ut: unittst

unittst: UnitTesting.o PTSymbolicObjects.o PathIntegration.o Multithreading.o Debugging.o FeynmanDiagram.o ExpressionSerialization.o
	$(CC) $(CFLAGS) UnitTesting.o PTSymbolicObjects.o PathIntegration.o Multithreading.o Debugging.o FeynmanDiagram.o ExpressionSerialization.o -o unittst $(LIBBOOST)
	
UnitTesting.o: UnitTesting.cpp
	$(CC) $(CFLAGS) -c UnitTesting.cpp

clean:
	rm *.o amaunet unittst
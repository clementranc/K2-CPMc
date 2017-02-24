#!/bin/bash
# ===================================================================#
# Makefile to compile libcpm, written for UNIX OS. It creates a shared
# library in C. The compilation commands can differ from one OS to
# another. This version should work if you use clang++ to compile C++.
# Some changes should be made to use GCC.
# ===================================================================#
#
# Filename of the targets
CT1=libcpm
CT2=table
CT3=matrix

# C++ compiler
CC=clang++

# C++ flags
CCFLAGS1:=-Wall -fPIC
CCFLAGS2:=-shared -undefined dynamic_lookup

# C++ libraries and linkers
#CCLIBS:=-L/opt/local/lib/gcc6/ -lgfortran
CCLINKERS:=-I/Applications/anaconda/include/python2.7 -I/Applications/anaconda/lib/python2.7

# Compilation instructions
libcpm:
	$(CC) $(CCFLAGS1) -c $(CT1).cpp $(CT2).cpp $(CT3).cpp
	swig -c++ -python $(CT1).i
	$(CC) $(CCFLAGS1) $(CCLINKERS) -c $(CT1)_wrap.cxx
	$(CC) $(CCFLAGS2) -o _$(CT1).so $(CT1).o $(CT1)_wrap.o $(CT2).o $(CT3).o

clean:
	rm $(CT1).py $(CT1).o $(CT1)_wrap.cxx $(CT1)_wrap.o

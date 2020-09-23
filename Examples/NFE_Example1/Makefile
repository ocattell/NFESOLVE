### Makefile for NFE_Example1 ###
PROGRAM_NAME := NFE_Example1

# Dependencies
DEPEND := NFE_Example1_Driver.cpp NeuralField.cpp

# Set path to NFESOLVE directory (Change to wherever you have it installed)
NFESOLVE_DIR := $(HOME)/Documents/NFESOLVE/

NFESOLVE_INCDIR := $(NFESOLVE_DIR)/include
NFESOLVE_LIBDIR := $(NFESOLVE_DIR)/lib

# Set up
CC := g++
CFLAGS := -std=c++11 -O2
OMP := -fopenmp
INC := -I $(NFESOLVE_INCDIR)
LIB_INCLUDES := -larmadillo


all: $(PROGRAM_NAME)

$(PROGRAM_NAME) : $(DEPEND)
	$(CC) $(CFLAGS) $(INC) -o $(PROGRAM_NAME) $^ -L$(NFESOLVE_LIBDIR) -lNFESOLVE $(LIB_INCLUDES)

clean :
	rm -f $(PROGRAM_NAME)

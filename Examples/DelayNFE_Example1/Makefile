### Makefile for DelayNFE_Example1 ###
PROGRAM_NAME := DelayNFE_Example1

# Dependencies
DEPEND := DelayNFE_Example1_Driver.cpp DelayNeuralField.cpp

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


$(PROGRAM_NAME) : $(DEPEND)
	$(CC) $(CFLAGS) $(INC) -o $(PROGRAM_NAME) $^ -L$(NFESOLVE_LIBDIR) -lNFESOLVE $(LIB_INCLUDES)

$(PROGRAM_NAME)_PAR : $(DEPEND)
	$(CC) $(CFLAGS) $(OMP) $(INC) -o $(PROGRAM_NAME) $^ -L$(NFESOLVE_LIBDIR) -lNFESOLVE_PAR $(LIB_INCLUDES)

clean :
	rm -f $(PROGRAM_NAME)

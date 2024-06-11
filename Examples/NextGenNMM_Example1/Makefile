### Makefile for NextGenNMM ###
PROGRAM_NAME := NextGenNMM

# Set path to NFESOLVE directory (Change to wherever you have it installed)
NFESOLVE_DIR := $(HOME)/NFESOLVE

NFESOLVE_INCDIR := $(NFESOLVE_DIR)/include
NFESOLVE_LIBDIR := $(NFESOLVE_DIR)/lib

# Dependencies
DEPEND_7 := NextGenNMM_Delay_Sparse_BOLD_Driver.cpp NextGenNMM_Delay_Sparse_BOLD.cpp

# Set up
CC := g++
CFLAGS := -std=c++11 -O2
OMP := -fopenmp
INC := -I $(NFESOLVE_INCDIR)
LIB_INCLUDES := -larmadillo


$(PROGRAM_NAME)_Delay_Sparse_BOLD : $(DEPEND_7)
	$(CC) $(CFLAGS) $(INC) -o $(PROGRAM_NAME).out $^ -L$(NFESOLVE_LIBDIR) -lNFESOLVE $(LIB_INCLUDES)

$(PROGRAM_NAME)_Delay_Sparse_BOLD_PAR : $(DEPEND_7)
	$(CC) $(CFLAGS) $(OMP) $(INC) -o $(PROGRAM_NAME).out $^ -L$(NFESOLVE_LIBDIR) -lNFESOLVE_PAR $(LIB_INCLUDES)

clean :
	rm -f $(PROGRAM_NAME).out

PPMDIR	= /home/kevin/ppm-library/ppm
PPMNUMDIR = /home/kevin/ppm-library/ppmnumerics
HDF5 = /home/kevin/Downloads/hdf5-1.8.11/hdf5

FC = mpif90
LDFLAGS = -L/home/kevin/ppm-library/metis/lib -latlas -lblas -lmetis -L$(PPMDIR)/lib -L$(PPMNUMDIR)/lib -lppmnumerics -lppm -lmetis -llapack -lblas -L$(HDF5)/lib -lhdf5_fortran -lhdf5hl_fortran

OPT= -O3 -ffree-form -ffree-line-length-none -ffree-line-length-none -O0 -g -fbounds-check -fbacktrace -ftrapv -fno-automatic -ffree-line-length-none -Wall -ffpe-trap=zero,overflow,underflow -I/home/kevin/ppm-library/checkpoint -I/usr/include -lhdf5_fortran

INC = -I$(PPMDIR)/include/ppm -I$(PPMNUMDIR)/include/ppm -I$(HDF5)/include
FFLAGS = $(OPT) $(INC) -cpp

SOURCES=pointer_tracker.f checkpoint_hdf5.f
OBJECTS=pointer_tracker.o checkpoint_hdf5.o

COMPILECMD = $(FC) $(INC) $(FFLAGS) -c -o $@ $<

.PHONY: source lib

all: test

test: lib test.f
	$(FC) -o test test.f $(INC) -L. -lcheck $(LDFLAGS)
	$(FC) -o polytest polytest.f $(INC) -L. -lcheck $(LDFLAGS)

source:
	$(MAKE) -C pointers
	$(MAKE) -C types
	$(MAKE) -C trees

lib: source clean $(OBJECTS)
	ar r libcheck.a $(OBJECTS)

clean:
	rm -f *.o

.f.o:
	$(COMPILECMD)

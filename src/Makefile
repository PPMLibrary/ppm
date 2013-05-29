PPMDIR	= /home/kevin/ppm-library/ppm
PPMNUMDIR = /home/kevin/ppm-library/ppmnumerics

FC = mpif90
LDFLAGS = -L/home/kevin/ppm-library/metis/lib -latlas -lblas -lmetis -L$(PPMDIR)/lib -L$(PPMNUMDIR)/lib -lppmnumerics -lppm -lmetis -llapack -lblas -lhdf5_fortran -lhdf5hl_fortran

OPT= -O3 -ffree-form -ffree-line-length-none -ffree-line-length-none -O0 -g -fbounds-check -fbacktrace -ftrapv -fno-automatic -ffree-line-length-none -Wall -ffpe-trap=zero,overflow,underflow -I/home/kevin/ppm-library/checkpoint -I/usr/include -lhdf5_fortran

INC = -I$(PPMDIR)/include/ppm -I$(PPMNUMDIR)/include/ppm -I/usr/include
FFLAGS = $(OPT) $(INC)

SOURCES=checkpoint_hdf5.f
OBJECTS=checkpoint_hdf5.o

COMPILECMD = $(FC) $(INC) $(FFLAGS) -c -o $@ $<

all: test

test: lib test.f
	$(FC) -o test test.f $(INC) -L. -lcheck $(LDFLAGS)
	

lib: clean $(OBJECTS)
	ar r libcheck.a $(OBJECTS)

clean:
	rm -f *.o

.f.o:
	$(COMPILECMD)

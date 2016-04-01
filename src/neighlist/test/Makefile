TESTINC=-I../../../include -L../../../lib -lppm
FC = mpif90 -ffree-form -ffree-line-length-none -ffree-line-length-none
sources=$(wildcard *.f)
targets=$(patsubst %.f,%,$(sources))
LDFLAGS = -L/Users/afshar/Documents/MPI_CBG/mosaic/code/metis4/lib -framework Accelerate
CPP = cpp-4.9

new : clean all

all : $(targets)

$(targets) : $(sources)
	rm -f __$@.f
	$(CPP) -P $@.f > __$@.f
	$(FC) $(OPTS) __$@.f $(TESTINC) $(LDFLAGS) $(LIBS) -o $@
	rm -f __$@.f
	./$@


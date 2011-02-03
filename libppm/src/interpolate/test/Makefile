TESTINC=-I../../../include -L../../../lib -lppm
sources=$(wildcard *.f)
targets=$(patsubst %.f,%,$(sources))

new : clean all run

all : $(targets)

$(targets) : $(sources)
	$(CPP) -P $@.f > __$@.f
	$(FC) $(OPTS) __$@.f $(TESTINC) $(LDFLAGS) $(LIBS) -o $@
	rm __$@.f
	./$@

clean :
	rm -f $(targets)


all: ppm_pointers.f

ppm_pointers.f: pointer.in typelist.in
	./gen_get_ptr.pl $@

clean:
	rm *.f

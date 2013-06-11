SOURCES = ppm_pointers.f ppm_t_ptr.f

all: $(SOURCES)

ppm_t_ptr.f: ptrlist.in ppm_t_ptr.in
	./store_ptr_list.pl

ppm_pointers.f: pointer.in typelist.in
	./gen_get_ptr.pl $@

clean:
	rm *.f

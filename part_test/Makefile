.PHONY: all

all:
	$(MAKE) -C gen

test: all
	mpirun -np 2 bin/part_test

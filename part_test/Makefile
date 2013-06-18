.PHONY: all

all:
	$(MAKE) -C gen

test: all
	mpirun -np 1 bin/part_test

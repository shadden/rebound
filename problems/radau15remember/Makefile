PROBLEMDIR=$(shell basename `dirname \`pwd\``)"/"$(shell basename `pwd`)
export OPENGL=1
export OPENMP=0
export OPT=-O3
export CC=gcc 
export QUADRUPOLE=0
export MPI=0
export LIB=-lgmp

radau15:
	# Setup link to different modules
	ln -fs gravity_direct.c ../../src/gravity.c
	ln -fs boundaries_open.c ../../src/boundaries.c
	ln -fs collisions_none.c ../../src/collisions.c
	ln -fs integrator_radau15.c ../../src/integrator.c
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/nbody .

leapfrog:
	# Setup link to different modules
	#ln -fs gravity_tree.c ../../src/gravity.c
	ln -fs gravity_direct.c ../../src/gravity.c
	ln -fs boundaries_open.c ../../src/boundaries.c
	ln -fs collisions_none.c ../../src/collisions.c
	ln -fs integrator_leapfrog.c ../../src/integrator.c
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/nbody .

doc: all
	$(MAKE) -C ../../src/ doc

clean:
	$(MAKE) -C ../../src/ clean
	rm -vf nbody

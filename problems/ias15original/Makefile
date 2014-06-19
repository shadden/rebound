PROBLEMDIR=$(shell basename `dirname \`pwd\``)"/"$(shell basename `pwd`)
export OPENGL=0
export OPT=-O3
export LIB=

ias15:
	# Setup link to different modules
	ln -fs gravity_direct.c ../../src/gravity.c
	ln -fs boundaries_open.c ../../src/boundaries.c
	ln -fs collisions_none.c ../../src/collisions.c
	ln -fs integrator_ias15.c ../../src/integrator.c
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/nbody .

ias15original:
	# Setup link to different modules
	#ln -fs gravity_tree.c ../../src/gravity.c
	ln -fs gravity_direct.c ../../src/gravity.c
	ln -fs boundaries_open.c ../../src/boundaries.c
	ln -fs collisions_none.c ../../src/collisions.c
	ln -fs ../$(PROBLEMDIR)/integrator_ias15original.c ../../src/integrator.c
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/nbody .


PROBLEMDIR=problems/$(shell basename `pwd`)
export OPT=-O3
export OPENGL=0
#export QUADRUPOLE=1
#export MPI=1
#export OPENMP=1
#export LIBPNG=1

# Note to myself:
#  - The L.H.S. of the ln commands set the *specific* code file to use
#  - The R.H.S. is a sort *variable* name that allows different modules to be 
#	used  
all:
	# Setup link to different modules
	ln -fs gravity_direct.c ../../src/gravity.c
	ln -fs boundaries_open.c ../../src/boundaries.c
	ln -fs integrator_ias15.c ../../src/integrator.c
	ln -fs collisions_none.c ../../src/collisions.c
	# Setup link to problem file
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

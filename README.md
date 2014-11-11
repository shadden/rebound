REBOUND - An open-source multi-purpose N-body code for collisional dynamics
===========================================================================

[![GPL](http://img.shields.io/badge/license-GPL-green.svg?style=flat)](https://github.com/hannorein/rebound/blob/master/LICENSE)
[![arXiv](http://img.shields.io/badge/arXiv-1110.4876-orange.svg?style=flat)](http://arxiv.org/abs/1110.4876)
[![arXiv](http://img.shields.io/badge/arXiv-1409.4779-orange.svg?style=flat)](http://arxiv.org/abs/1409.4779)


Contributors
------------
* Hanno Rein, University of Toronto, <hanno@hanno-rein.de>
* Shangfei Liu, Kavli Institute for Astronomy and Astrophysics at Peking University (KIAA-PKU), Beijing, <liushangfei@pku.edu.cn>
* David S. Spiegel, Institute for Advanced Study (IAS), Princeton, <dave@ias.edu>
* Akihiko Fujii, National Astronomical Observatory of Japan/University of Tokyo, Tokyo, <akihiko.fujii@nao.ac.jp>

Papers
----- 

There are two papers describing the functionality of REBOUND. 

The first one, Rein & Liu (Astronomy and Astrophysics, Volume 537, A128, 2012, http://arxiv.org/abs/1110.4876), describes the code structure and the main feature including the gravity and collision routines. 

The second paper, Rein & Spiegel (MNRAS, in press, http://arxiv.org/abs/1409.4779) describes the versatile high order integrator IAS15 which is now part of REBOUND. 


Screenshot
---------- 
 
![Tree structure in REBOUND](https://raw.github.com/hannorein/rebound/master/doc/images/screenshot_shearingsheet.png) 

You can also find a video on YouTube, http://youtu.be/gaExPGW1WzI?hd=1, that shows how to download and install REBOUND. 


Python and other programming languages 
-----------------
REBOUND is written in C. However, we provide a simple dynamic library `libias15` for the new IAS15 integrator. Most of the features that make REBOUND great are not available in this library. It is purely for efficiently accessing the high order integrator IAS15. There is no visualization and no helper function to setup particles. To compile this library, go to the `shared` folder and type `make`.

The most interesting use case for `libias15` is a python wrapper that we provide. This wrapper can be used to very easily access `libias15`. The wrapper (module) might appeal to people who want to setup their problem in python and then call IAS15 to efficiently integrate particles with very high precision. The following listing shows a complete python script to run an N-body simulation with IAS15 and `libias15`:
 
```python
# Import the rebound module
import sys; sys.path.append('../')
import rebound

# Add particles
rebound.particle_add( rebound.Particle(m=1.) )                  # Star
rebound.particle_add( rebound.Particle(m=1e-3,x=1.,vy=1.) )     # Planet

# Move particles so that the center of mass is (and stays) at the origin  
rebound.move_to_center_of_momentum()

# Integrate until t=100 (roughly 16 orbits) 
rebound.integrate(100.)
```

For details, please look at the README file in the `python_examples` directory. 


Available modules
-----------------

REBOUND is extremely modular. You have the choice between different gravity, collision, boundary and integration modules. It is also possible to implement completely new modules with minimal effort. Modules are chosen by setting symbolic links. Thus, there is no need to run a configure script. For example, there is one link `gravity.c` that points to one of the gravity modules `gravity_*.c`. The symbolic links are set in the problem makefile (see below).

This setup allows you to work on multiple projects at the same time using different modules. When switching to another problem, nothing has to be set-up and the problem can by compiled by simply typing `make` in the corresponding directory (see below).

### Gravity ###
<table>
  <tr><th>Module name</th>
     <th>Description</th></tr>
  <tr><td><pre>gravity_none.c       </pre></td>
     <td>No self-gravity</td></tr>
  <tr><td><pre>gravity_direct.c     </pre></td>
     <td>Direct summation, O(N^2)</td></tr>
  <tr><td><pre>gravity_opencl.c     </pre></td>
     <td>Direct summation, O(N^2), but accelerated using the OpenCL framework. </td></tr>
  <tr><td><pre>gravity_tree.c       </pre></td>
     <td>Oct tree, Barnes & Hut 1986, O(N log(N))</td></tr>
  <tr><td><pre>gravity_grape.c      </pre></td>
     <td>GRAPE, hardware accelerated direct summation, Sugimoto et al. 1990 </td></tr>
  <tr><td><pre>gravity_fft.c        </pre></td>
     <td>Two dimensional gravity solver using FFTW, works in a periodic box and the shearing sheet. (Not well tested yet.)</td></tr>
</table>

### Collision detection ###
<table>
  <tr><th>Module name</th>
     <th>Description</th></tr>
  <tr><td><pre>collisions_none.c    </pre></td>
     <td>No collision detection</td></tr>
  <tr><td><pre>collisions_direct.c  </pre></td>
     <td>Direct nearest neighbor search, O(N^2)</td></tr>
  <tr><td><pre>collisions_tree.c    </pre></td>
     <td>Oct tree, O(N log(N))</td></tr>
  <tr><td><pre>collisions_sweep.c   </pre></td>
     <td>Plane sweep algorithm, ideal for low dimensional problems, O(N) or O(N^1.5) depending on geometry</td></tr>
  <tr><td><pre>collisions_sweepphi.c   </pre></td>
     <td>Plane sweep algorithm along the azimuthal angle, ideal for narrow rings in global simulations, O(N) or O(N^1.5) depending on geometry</td></tr>
</table>

### Integrators ###
<table>
  <tr><th>Module name</th>
     <th>Description</th></tr>
  <tr><td><pre>integrator_euler.c   </pre></td>
     <td>Euler scheme, first order</td></tr>
  <tr><td><pre>integrator_leapfrog.c</pre></td>
     <td>Leap frog, second order, symplectic</td></tr>
  <tr><td><pre>integrator_wh.c      </pre></td>
     <td>Wisdom-Holman Mapping, mixed variable symplectic integrator for the Kepler potential, second order, Wisdom & Holman 1991, Kinoshita et al 1991</td></tr>
  <tr><td><pre>integrator_ias15.c </pre></td>
     <td>IAS15 stands for Integrator with Adaptive Step-size control, 15th order. It is a vey high order, non-symplectic integrator which can handle arbitrary (velocity dependent) forces and is in most cases accurate down to machine precission. Rein & Spiegel 2014, Everhart 1985</td></tr>
  <tr><td><pre>integrator_sei.c     </pre></td>
     <td>Symplectic Epicycle Integrator (SEI), mixed variable symplectic integrator for the shearing sheet, second order, Rein & Tremaine 2011</td></tr>
</table>

### Boundaries ###
<table>
  <tr><th>Module name</th>
     <th>Description</th></tr>
  <tr><td><pre>boundaries_open.c    </pre></td>
     <td>Particles are removed from the simulation if they leaves the box.</td></tr>
  <tr><td><pre>boundaries_periodic.c</pre></td>
     <td>Periodic boundary conditions. Particles are reinserted on the other side if they cross the box boundaries. You can use an arbitrary number of ghost-boxes with this module.</td></tr>
  <tr><td><pre>boundaries_shear.c   </pre></td>
     <td>Shear periodic boundary conditions. Similar to periodic boundary conditions, but ghost-boxes are moving with constant speed, set by the shear.</td></tr>
</table>


Other features worth mentioning
-------------------------------
* Real-time, 3D OpenGL visualization.
* The code is written entirely in C. It conforms to the ISO standard C99.
* Parallelized with OpenMP (for shared memory systems).
* Parallelized with MPI using an essential tree for gravity and collisions (for distributed memory systems).
* No libraries are needed. The use of OpenGL/GLUT/libpng for visualization is optional. 
* The code is fully open-source and can be downloaded freely from http://github.com/hannorein/rebound.
* No configuration is needed to run any of the example problems. Just type `make && ./nbody` in the problem directory to run them.
* Standard ASCII or binary output routines. 
* Different modules are easily interchangeable by one line in the Makefile.
  

How to download, compile and run REBOUND
----------------------------------------

### For the impatient ###
If you are using a Mac, make sure you have a compiler suite installed. Open a terminal and type `make`. If it is not installed, go to the AppStore and download Xcode (it is free). Once installed, open Xcode, go to Settings, then Downloads and install the Command Line Tools. 

Then, simply copy and paste this line to your terminal and press enter

    git clone http://github.com/hannorein/rebound && cd rebound/examples/shearing_sheet && make && ./nbody

or if you do not have git installed

    wget --no-check-certificate https://github.com/hannorein/rebound/tarball/master -O- | tar xvz && cd hannorein-rebound-*/examples/shearing_sheet/ && make && ./nbody

### For the patient ###
REBOUND is very easy to install and use. To get started, download the latest version of the code from github. If you are familiar with `git`, you can clone the project and keep up-to-date with the latest developments. Otherwise, you can also simply download a snapshot of the repository as a tar or zip file at http://github.com/hannorein/rebound. There is a download bottom at the top right. 

In the main directory, you find a sub-directory called `src` which contains the bulk parts of the  source code and a directory called `examples` with various example problems. To compile one of the example, you have to go to that directory, for example:

    cd examples/shearing_sheet/

Then, type

    make

This will do the following things    

* It sets various environment variables. These determine settings like the compiler optimization flags and which libraries are included (see below). 
* It creates symbolic links to the active modules. This allows you to choose from different gravity solvers, boundary conditions, integrators and collision solvers. For example, to change the gravity solver from using a tree to direct summation you could change `gravity_tree.c` to `gravity_direct.c`. 
* It creates a symbolic link to the current problem file. Each problem file contains the initial conditions and the output routines for the current problem. You do not need to change any file in `src/` to create a new problem unless you want to do something very special. This keeps the initial conditions and the code itself cleanly separated.
* It compiles the code and copies the binary into the current directory.

If something goes wrong, it is most likely the visualization module. You can turn it off by deleting the line which contains `OPENGL` in the makefile. Of course, you will not see much unless you put in some extra work to visualize the results.

You can also create a documentation based on the current choice of modules by typing `make doc`. However, this requires the documentation generator `doxygen` to be installed. The documentation will be generated in the directory `doc/html/`.

To finally run the code, simply type

    ./nbody

A window should open and you will see a simulation running in real time. The problem in the directory `examples/shearing_sheet/` simulates the rings of Saturn and uses a local shearing sheet approximation. Have a look at the other examples as well and you will quickly get an idea of what REBOUND can do. 

If you want to create your own problem, just copy one of the example directories or the template in the `problems` directory. Then simply modify `problem.c` and `Makefile` accordingly.  

### How to install GLUT ###
The OpenGL Utility Toolkit (GLUT) comes pre-installed as a framework on Mac OSX. If you are working on another operating system, you might have to install GLUT yourself if you see an error message such as `error: GL/glut.h: No such file or directory`. On Debian and Ubuntu, simply make sure the `freeglut3-dev` package is installed. If glut is not available in your package manager, go to http://freeglut.sourceforge.net/ download the latest version, configure it with `./configure` and compile it with `make`. Finally install the library and header files with `make install`. 

You can also install freeglut in a non-default installation directory if you do not have super-user rights by running the freeglut installation script with the prefix option:

    mkdir ${HOME}/local
    ./configure --prefix=${HOME}/local
    make all && make install

Then, add the following lines to the REBOUND Makefile

    OPT += -I$(HOME)/local/include
    LIB += -L$(HOME)/local/lib

Note that you can still compile and run REBOUND even if you do not have GLUT installed. Simple set `OPENGL=0` in the makefile (see below). 

### Environment variables ###
The makefile in each problem directory sets various environment variables. These determine the compiler optimization flags, the libraries included and basic code settings. Let us look at one of the examples `shearing_sheet` in more detail. 

- `export PROFILING=1`. This enables profiling. You can see how much time is spend in the collision, gravity, integrator and visualization modules. This is useful to get an idea about the computational bottleneck.
- `export QUADRUPOLE=0`. This disables the calculation of quadrupole moments for each cell in the tree. The simulation is faster, but less accurate.
- `export OPENGL=1`. This enables real-time OpenGL visualizations and requires both OpenGL and GLUT libraries to be installed. This should work without any further adjustments on any Mac which has Xcode installed. On Linux both libraries must be installed in `/usr/local/`. You can change the default search paths for libraries in the file `src/Makefile`. 
- `export MPI=0`. This disables parallelization with MPI.
- `export OPENMP=1`. This enables parallelization with OpenMP. The number of threads can be set with an environment variable at runtime, e.g.: `export OMP_NUM_THREADS=8`.
- `export CC=icc`. This flag can be used to override the default compiler. The default compilers are `gcc` for the sequential and `mpicc` for the parallel version. 
- `export LIB=`. Additional search paths for external libraries (such as OpenGL, GLUT and LIBPNG) can be set up using this variable. 
- `export OPT=-O3`. This sets the additional compiler flag `-O3` and optimizes the code for speed. Additional search paths to header files for external libraries (such as OpenGL, GLUT and LIBPNG) can be set up using this variable. 

All of these variables are read by the main makefile in the `src/` directory. The `OPENGL` variable, for example, is used to determine if the OpenGL and GLUT libraries should be included. If the variable is `1` the makefile also sets a pre-compiler macro with `-DOPENGL`. Note that because OPENGL is incompatible with MPI, when MPI is turned on (set to 1), OPENGL is automatically turned off (set to 0) in the main makefile.


### User-defined functions in the problem.c file ###
The problem.c file contains at least four functions. You do not need to implement all of them, a dummy might be enough. 

#### void problem_init(int argc, char* argv[]) ####
This routine is where you read command line arguments and set up your initial conditions. REBOUND does not come with a built-in functionality to read configuration files at run-time. You should see this as a feature. In REBOUND, you have one `problem.c` file for each problem. Thus, everything can be set within this file. There are, of course, situation in which you want to do something like a parameter space survey. In almost all cases, you vary a few parameters but rarely more, say 5. You can easily read these parameters from the command line.
 
Here is one example that reads in the first argument given to rebound as the box-size and sets a default value when no value is given:

```c
if (argc>1){
	boxsize = atof(argv[1]);
}else{
	boxsize = 100;
}
```

If you are still convinced that you need a configuration file, you are welcome to implement it yourself. This function is where you want to do that.    

#### void problem_additional_forces() ####
This is a function pointer which is called one or more times per time-step whenever the forces are updated. This is where you can implement all kind of things such as additional forces onto particles. 

The following lines of code implement a simple velocity dependent force.  `integrator_ias15.c` is best suited for this (see `examples/dragforce`):

```c
void velocity_dependent_force(){
	for (int i=1;i<N;i++){
		particles[i].ax -= 0.0000001 * particles[i].vx;
		particles[i].ay -= 0.0000001 * particles[i].vy;
		particles[i].az -= 0.0000001 * particles[i].vz;
	}
}
```

Make sure you set the function pointer in the `problem_init()` routine:

```c
	problem_additional_forces = velocity_dependent_force;
```

By default, all integrators assume that the forces are velocity dependent. If all forces acting on particles only depend on positions, you can set the following variable (defined in `integrator.h`) to `0` to speed up the calculation:

```c
	integrator_force_is_velocitydependent = 0;
```


#### void problem_output() ####
This function is called at the beginning of the simulation and at the end of each time-step. You can implement your output routines here. Many basic output functions are already implemented in REBOUND. See `output.h` for more details. The function `output_check(odt)` can be used to easily check if an output is needed after a regular interval. For example, the following code snippet outputs some timing statistics to the console every 10 time-steps:

```c
if (output_check(10.*dt)){
	output_timing();
}
```    
 
#### void problem_finish() ####
This function is called at the end of the simulation, when t >= tmax. This is the last chance to output any quantities before the program ends.


Examples
-----------------------
The following examples can all be found in the `examples` directory. 
Whatever you plan to do with REBOUND, chances are there is already an example available which you can use as a starting point.


*  **examples/bouncing_balls**

  This example is using the following modules:  
  `gravity_direct.c`
  `boundaries_periodic.c`
  `integrator_leapfrog.c`
  `collisions_direct.c`

  This example is a simple test of collision detection
  methods. To change the collision detection algorithm, you can replace
  the module collisions_direct.c to either collisions_tree.c or
  collisions_sweep.c in the Makefile.
  

*  **examples/bouncing_balls_corners**

  This example is using the following modules:  
  `gravity_direct.c`
  `boundaries_periodic.c`
  `integrator_leapfrog.c`
  `collisions_tree.c`

  This example tests collision detection methods accros box boundaries.
  There are four particles, one in each corner. To see the ghost boxes in OpenGL
  press `g` while the simulation is running.
  

*  **examples/bouncing_string**

  This example is using the following modules:  
  `gravity_none.c`
  `boundaries_periodic.c`
  `integrator_leapfrog.c`
  `collisions_direct.c`

  This example tests collision detection methods.
  The example uses a non-square, rectangular box. 10 particles are placed
  along a line. All except one of the particles are at rest
  initially.
  

*  **examples/circumplanetarydust**

  This example is using the following modules:  
  `gravity_direct.c`
  `boundaries_open.c`
  `integrator_ias15.c`
  `collisions_none.c`

  This example shows how to integrate circumplanetary
  dust particles using the `integrator_ias15.c` module.
  The example sets the function pointer `problem_additional_forces`
  to its own function that describes the radiation forces.
  The example uses a beta parameter of 0.01.
  The output is custom too, outputting the semi-major axis of
  every dust particle relative to the planet.
  Only one dust particle is used in this example, but there could be
  many.
  

*  **examples/closeencounter**

  This example is using the following modules:  
  `gravity_direct.c`
  `boundaries_open.c`
  `integrator_ias15.c`
  `collisions_none.c`

  This example integrates a densly packed planetary system
  which becomes unstable on a timescale of only a few orbits. The IAS15
  integrator with adaptive timestepping is used. This integrator
  automatically decreases the timestep whenever a close
  enocunter happens. IAS15 is very high order and ideally suited for the
  detection of these kind of encounters.
  

*  **examples/closeencounter_record**

  This example is using the following modules:  
  `gravity_direct.c`
  `boundaries_open.c`
  `integrator_ias15.c`
  `collisions_direct.c`

  This example integrates a densly packed planetary system
  which becomes unstable on a timescale of only a few orbits.
  The example is identical to the `close_encounter` sample, except that
  the collisions are recorded and written to a file. What kind of collisions
  are recorded can be easily modified. It is also possible to implement some
  additional physics whenever a collision has been detection (e.g. fragmentation).
  The collision search is by default a direct search, i.e. O(N^2) but can be
  changed to a tree by using the `collisions_tree.c` module.
  

*  **examples/dragforce**

  This example is using the following modules:  
  `gravity_none.c`
  `boundaries_open.c`
  `collisions_none.c`
  `integrator_ias15.c`

  This is a very simple example on how to implement a velocity
  dependent drag force. The example uses the IAS15 integrator, which
  is ideally suited to handle non-conservative forces.
  No gravitational forces or collisions are present.
  

*  **examples/eccentric_orbit**

  This example is using the following modules:  
  `gravity_direct.c`
  `boundaries_open.c`
  `integrator_ias15.c`
  `collisions_none.c`

  This example uses the IAS15 integrator to simulate
  a very eccentric planetary orbit. The integrator
  automatically adjusts the timestep so that the pericentre passages
  resovled with high accuracy.
  
  

*  **examples/forced_migration**

  This example is using the following modules:  
  `gravity_direct.c`
  `boundaries_open.c`
  `integrator_ias15.c`
  `collisions_none.c`

  This example applies dissipative forces to two
  bodies orbiting a central object. The forces are specified
  in terms of damping timescales for the semi-major axis and
  eccentricity. This mimics planetary micration in a protostellar disc.
  The example reproduces the study of Lee & Peale (2002) on the
  formation of the planetary system GJ876. For a comparison,
  see figure 4 in their paper. The IAS15 integrator is used
  because the forces are velocity dependent.
  Special thanks goes to Willy Kley for helping me to implement
  the damping terms as actual forces.
  
  

*  **examples/granulardynamics**

  This example is using the following modules:  
  `gravity_none.c`
  `boundaries_periodic.c`
  `integrator_leapfrog.c`
  `collisions_tree.c`

  This example is about granular dynamics. No gravitational
  forces are present in this example, which is why the module
  `gravity_none.c` is used. Two boundary layers made of
  particles simulate shearing walls. These walls are heating
  up the particles, create a dense and cool layer in the middle.
  

*  **examples/J2**

  This example is using the following modules:  
  `gravity_direct.c`
  `boundaries_open.c`
  `integrator_ias15.c`
  `collisions_none.c`

  This example presents an implementation of the J2
  gravitational moment. The equation of motions are integrated with
  the 15th order IAS15 integrator. The parameters in this examples
  have been chosen to represent those of Saturn, but you can easily
  change them or even include higher order terms in the multipole
  expansion.
  

*  **examples/kozai**

  This example is using the following modules:  
  `gravity_direct.c`
  `boundaries_open.c`
  `integrator_ias15.c`
  `collisions_none.c`

  This example uses the IAS15 integrator to simulate
  a Lidov Kozai cycle of a planet perturbed by a distant star. The integrator
  automatically adjusts the timestep so that even very high
  eccentricity encounters are resovled with high accuracy.
  
  

*  **examples/mergers**

  This example is using the following modules:  
  `gravity_direct.c`
  `boundaries_open.c`
  `integrator_ias15.c`
  `collisions_direct.c`

  This example integrates a densly packed planetary system
  which becomes unstable on a timescale of only a few orbits. The IAS15
  integrator with adaptive timestepping is used. The bodies have a finite
  size and merge if they collide. Note that the size is unphysically large
  in this example.
  

*  **examples/opencl**

  This example is using the following modules:  
  `gravity_opencl.c`
  `boundaries_open.c`
  `integrator_leapfrog.c`
  `collisions_none.c`
  `gravity_direct.c`
  `boundaries_open.c`
  `integrator_leapfrog.c`
  `collisions_none.c`

  A self-gravitating disc is integrated using
  the OpenCL direct gravity summation module.
  
  This is a very simple implementation (see `gravity_opencl.c`).
  Currently it only supports floating point precission. It also
  transfers the data back and forth from the GPU every timestep.
  There are considerable improvements to be made. This is just a
  proof of concept. Also note that the code required N to be a
  multiple of the workgrop size.
  
  You can test the performance increase by running:
  `make direct && ./nbody`, which will run on the CPU and
  `make && ./nbody`, which will run on the GPU.
  
  The Makefile is working with the Apple LLVM compiler. Changes
  might be necessary for other compilers such as gcc.
  
  

*  **examples/outer_solar_system**

  This example is using the following modules:  
  `gravity_direct.c`
  `boundaries_open.c`
  `integrator_ias15.c`
  `collisions_none.c`

  This example uses the IAS15 integrator
  to integrate the outer planets of the solar system. The initial
  conditions are taken from Applegate et al 1986. Pluto is a test
  particle. This example is a good starting point for any long term orbit
  integrations.
  
  You probably want to turn off the visualization for any serious runs.
  Just go to the makefile and set `OPENGL=0`.
  
  The example also works with the Wisdom-Holman symplectic integrator.
  Simply change the integrator to `integrator_wh.c` in the Makefile.
  

*  **examples/overstability**

  This example is using the following modules:  
  `gravity_none.c`
  `boundaries_shear.c`
  `integrator_sei.c`
  `collisions_sweep.c`

  A narrow box of Saturn's rings is simulated to
  study the viscous overstability. Collisions are resolved using
  the plane-sweep method.
  
  It takes about 30 orbits for the overstability to occur. You can
  speed up the calculation by turning off the visualization. Just press
  `d` while the simulation is running. Press `d` again to turn it back on.
  
  You can change the viewing angle of the camera with your mouse or by pressing
  the `r` key.
  

*  **examples/prdrag**

  This example is using the following modules:  
  `gravity_direct.c`
  `boundaries_open.c`
  `integrator_ias15.c`
  `collisions_none.c`

  This example provides an implementation of the
  Poynting-Robertson effect. The code is using the IAS15 integrator
  which is ideally suited for this velocity dependent force.
  

*  **examples/restarting_simulation**

  This example is using the following modules:  
  `gravity_direct.c`
  `boundaries_shear.c`
  `integrator_sei.c`
  `collisions_direct.c`

  This example demonstrates how to restart a simulation
  using a binary file. A shearing sheet ring simulation is used, but
  the same method can be applied to any other type of simulation.
  
  First, run the program with `./nbody`.
  Random initial conditions are created and
  a restart file is written once per orbit.
  Then, to restart the simulation, run the
  program with `./nbody --restart restart.bin`.
  
  

*  **examples/restricted_threebody**

  This example is using the following modules:  
  `gravity_direct.c`
  `boundaries_open.c`
  `integrator_wh.c`
  `collisions_none.c`

  This example simulates a disk of test particles around
  a central object, being perturbed by a planet.
  

*  **examples/restricted_threebody_mpi**

  This example is using the following modules:  
  `gravity_direct.c`
  `boundaries_open.c`
  `integrator_leapfrog.c`
  `collisions_none.c`

  This problem uses MPI to calculate the restricted three
  body problem. Active particles are copied to all nodes. All other
  particles only exist on one node and are not automatically (re-)
  distributed. There is not domain decomposition used in this example.
  Run with `mpirun -np 4 nbody`.
  

*  **examples/selfgravity_disc**

  This example is using the following modules:  
  `gravity_tree.c`
  `boundaries_open.c`
  `integrator_leapfrog.c`
  `collisions_none.c`

  A self-gravitating disc is integrated using
  the leap frog integrator. This example is also compatible with
  the Wisdom Holman integrator or the IAS15 integrator. Collisions
  are not resolved.
  

*  **examples/selfgravity_disc_grape**

  This example is using the following modules:  
  `gravity_grape.c`
  `boundaries_open.c`
  `integrator_leapfrog.c`
  `collisions_none.c`

  A self-gravitating disc is integrated using
  the leap frog integrator. This example is using the GRAPE
  module to calculate the self-gravity. You need to have a physical
  GRAPE card in your computer to run this example.
  Collisions are not resolved.
  

*  **examples/selfgravity_plummer**

  This example is using the following modules:  
  `gravity_tree.c`
  `boundaries_open.c`
  `integrator_leapfrog.c`
  `collisions_none.c`

  A self-gravitating plummer sphere is integrated using
  the leap frog integrator. Collisions are not resolved. Note that the
  fixed timestep might not allow you to resolve individual two-body
  encounters. An alternative integrator is `integrator_ias15.c` which
  comes with adaptive timestepping.
  

*  **examples/shearing_sheet**

  This example is using the following modules:  
  `gravity_tree.c`
  `boundaries_shear.c`
  `integrator_sei.c`
  `collisions_tree.c`

  This example simulates a small patch of Saturn's
  Rings in shearing sheet coordinates. If you have OpenGL enabled,
  you'll see one copy of the computational domain. Press `g` to see
  the ghost boxes which are used to calculate gravity and collisions.
  Particle properties resemble those found in Saturn's rings.
  
  

*  **examples/shearing_sheet_2**

  This example is using the following modules:  
  `gravity_tree.c`
  `boundaries_shear.c`
  `integrator_sei.c`
  `collisions_tree.c`

  This example is identical to the shearing_sheet
  example but uses a different algorithm for resolving individual
  collisions. In some cases, this might give more realistic results.
  Particle properties resemble those found in Saturn's rings.
  
  In this collision resolve method, particles are displaced if they
  overlap. This example also shows how to implement your own collision
  routine. This is where one could add fragmentation, or merging of
  particles.
  
  

*  **examples/shearing_sheet_fft**

  This example is using the following modules:  
  `gravity_fft.c`
  `boundaries_shear.c`
  `integrator_sei.c`
  `collisions_sweep.c`

  This problem is identical to the other shearing
  sheet examples but uses an FFT based gravity solver.
  To run this example, you need to install the FFTW library.
  Collisions are detected using a plane sweep algorithm.
  There is no tree present in this simulation.
  

*  **examples/shearing_sheet_grape**

  This example is using the following modules:  
  `gravity_grape.c`
  `boundaries_shear.c`
  `integrator_sei.c`
  `collisions_sweep.c`

  This is yet another shearing sheet example,
  it uses a GRAPE to calculate gravity. Note that you need to have
  a physical GRAPE card installed in your computer to run this
  simulation. Particle properties resemble those found in
  Saturn's rings.
  

*  **examples/shearing_sheet_profiling**

  This example is using the following modules:  
  `gravity_tree.c`
  `boundaries_shear.c`
  `integrator_sei.c`
  `collisions_tree.c`

  This example demonstrates how to use the
  profiling tool that comes with REBOUND to find out which parts
  of your code are slow. To turn on this option, simple set
  `PROFILING=1` in the Makefile.
  

*  **examples/spreading_ring**

  This example is using the following modules:  
  `gravity_direct.c`
  `boundaries_open.c`
  `integrator_wh.c`
  `collisions_sweepphi.c`

  A narrow ring of collisional particles is spreading.
  The example uses the Wisdom Holman integrator. A plane-sweep algorithm
  in the phi direction is used to detect collisions.
  

*  **examples/symplectic_integrator**

  This example is using the following modules:  
  `gravity_direct.c`
  `boundaries_open.c`
  `integrator_wh.c`
  `collisions_none.c`

  This example uses the symplectic Wisdom Holman (WH) integrator
  to integrate test particles on eccentric orbits in a fixed potential.
  Note that the WH integrator assumes that the central object is at the origin.
  

*  **examples/viewer**

  This example is using the following modules:  
  `gravity_none.c`
  `boundaries_periodic.c`
  `integrator_dummy.c`
  `collisions_dummy.c`

  This example doesn't simulate anything. It's just a
  visualization toll that can display data in the form x, y, z, r.
  This might be useful when large simulations have been run and you want
  to look (at parts of) it at a later time.
  
  Note that this example uses only dummy modules.
  
This viewer can display data in the form x, y, z, r.  
OpenGL keyboard command
-----------------------
You can use the following keyboard command to alter the OpenGL real-time visualizations.
<table>
  <tr><th>Key</th>
     <th>Description</th></tr>
  <tr><td>(space)</td><td>Pause simulation.</td></tr>
  <tr><td>d</td><td>Pause real-time visualization (simulation continues).</td></tr>
  <tr><td>q</td><td>Quit simulation.</td></tr>
  <tr><td>s</td><td>Toggle three dimensional spheres (looks better)/points (draws faster)</td></tr>
  <tr><td>g</td><td>Toggle ghost boxes</td></tr>
  <tr><td>r</td><td>Reset view. Press multiple times to change orientation.</td></tr>
  <tr><td>x/X</td><td>Move to a coordinate system centered on a particle (note: does not work if particle array is constantly resorted, i.e. in a tree.)</td></tr>
  <tr><td>t</td><td>Show tree structure.</td></tr>
  <tr><td>m</td><td>Show center of mass in tree structure (only available when t is toggled on).</td></tr>
  <tr><td>p</td><td>Save screen shot to file.</td></tr>
  <tr><td>c</td><td>Toggle clear screen after each time-step.</td></tr>
  <tr><td>w</td><td>Draw orbits as wires (particle with index 0 is central object).</td></tr>
</table>

Support and contributions
-------------------------
We offer limited support for REBOUND. If you encounter any problems, just send us an e-mail with as much details as possible and include your problem.c and makefile. Please make sure you are using the latest version of REBOUND that is available on github. 

REBOUND is open source and you are strongly encouraged to contribute to this project if you are using it. Please contact us and we will give you permission to push directly to the public repository. 


License
-------
REBOUND is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

REBOUND is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with REBOUND.  If not, see <http://www.gnu.org/licenses/>.


Acknowledgments
---------------
When you use this code or parts of this code for results presented in a scientific publication, please send us a copy of your paper so that we can keep track of all publications that made use of the code. We would greatly appreciate a citation to Rein and Liu (2012) and an acknowledgment of the form: 

_Simulations in this paper made use of the collisional N-body code REBOUND which can be downloaded freely at http://github.com/hannorein/rebound._

If you use the IAS15 integrator, please cite Rein and Spiegel (2014).

References in BibTeX format:

    @ARTICLE{ReinLiu2012,
       author = {{Rein}, H. and {Liu}, S.-F.},
        title = "{REBOUND: an open-source multi-purpose N-body code for collisional dynamics}",
      journal = {\aap},
    archivePrefix = "arXiv",
       eprint = {1110.4876},
     primaryClass = "astro-ph.EP",
     keywords = {methods: numerical, planets and satellites: rings, protoplanetary disks},
         year = 2012,
        month = jan,
       volume = 537,
          eid = {A128},
        pages = {A128},
          doi = {10.1051/0004-6361/201118085},
       adsurl = {http://adsabs.harvard.edu/abs/2012A%26A...537A.128R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

    @ARTICLE{ReinSpiegel2014,
       author = {{Rein}, H. and {Spiegel}, D.~S.},
        title = "{IAS15: A fast, adaptive, high-order integrator for gravitational dynamics, accurate to machine precision over a billion orbits}",
      journal = {ArXiv e-prints},
    archivePrefix = "arXiv",
       eprint = {1409.4779},
     primaryClass = "astro-ph.EP",
     keywords = {Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - Solar and Stellar Astrophysics, Mathematics - Numerical Analysis},
         year = 2014,
        month = sep,
       adsurl = {http://adsabs.harvard.edu/abs/2014arXiv1409.4779R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }



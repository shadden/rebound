#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"
#include "communication_mpi.h"
#include "tree.h"
#include "tools.h"

extern double OMEGA;
extern double OMEGAZ;
extern double minimum_collision_velocity;

extern double (*coefficient_of_restitution_for_velocity)(double); 
double coefficient_of_restitution_bridges(double v); 

double moonlet_radius;

void problem_init(int argc, char* argv[]){
	// Setup constants
	int oversample = 4;
	root_nx = oversample; root_ny = 4*oversample; root_nz = 1;
	nghostx = 2; nghosty = 3; nghostz = 0; 		// Ghost boxes, the more, the better, the slower.
	N_active	= 1;				// Only moonlet has gravity
	N_tree_fixed	= 1;				// Moonlet keeps the index 0
	
	OMEGA 		= 0.00013143527;		// 1/s
	OMEGAZ		= OMEGA*3.6;			// Enhancecd vertical epicyclic frequency
	G 		= 6.67428e-11;			// N m^2 / kg^2
	softening 	= 1;				// m (not really needed)
	dt 		= 1e-3*2.*M_PI/OMEGA;		// s (s small fraction of the orbital time)
	tmax		= 1e2 *2.*M_PI/OMEGA;		// Stop simulation after 100 orbits


	double surfacedensity 		= 400; 		// kg / m^2
	double particle_density		= 400;		// kg / m^3
	double particle_radius_min 	= 2.5;		// m
	double particle_radius_max 	= 3.5;		// m
	double particle_radius_slope 	= -3;	
	if (argc>1){					
		moonlet_radius 		= atof(argv[1]);// Try to read moonletsize from command line (in m)
	}else{
		moonlet_radius 		= 25;		// otherwise use default
	}
	boxsize = 1./(double)oversample*pow(moonlet_radius,2);		// boxsize scales as the third power of the moonlet
	
	// Use Bridges et al. coefficient of restitution.
	coefficient_of_restitution_for_velocity = coefficient_of_restitution_bridges;
	minimum_collision_velocity = particle_radius_min*OMEGA*0.001;  // small fraction of the shear over a particle_radius to avoid numerical difficulties
	
	// Setup particle structures
	init_box();
	// Initial conditions
	// Moonlet
	struct particle moonlet;
	moonlet.x 	= 0; moonlet.y 		= 0; moonlet.z 		= 0;
	moonlet.vx 	= 0; moonlet.vy 	= 0; moonlet.vz 	= 0;
	moonlet.ax 	= 0; moonlet.ay		= 0; moonlet.az		= 0;
	moonlet.r 	= moonlet_radius;		
	moonlet.m 	= particle_density*4./3.*M_PI*moonlet_radius*moonlet_radius*moonlet_radius;
	particles_add(moonlet);

	double total_mass = surfacedensity*boxsize*boxsize;
	for(int i=0;i<root_n;i++){
#ifdef MPI
		if (communication_mpi_rootbox_is_local(i)==0) continue;
#endif
		// Ring particles
		double mass = 0;
		while(mass<total_mass){
			struct particle pt;
			int ri = i%root_nx;
			int rj = ((i-ri)/root_nx)%root_ny;
			int xmin = -boxsize_x/2.+(double)(ri)*boxsize;
			int ymin = -boxsize_y/2.+(double)(rj)*boxsize;
			pt.x 		= tools_uniform(xmin,xmin+boxsize);
			pt.y 		= tools_uniform(ymin,ymin+boxsize);
			pt.z 		= tools_normal(1.);					// m
			if (pt.x*pt.x + pt.y*pt.y + pt.z*pt.z <moonlet_radius*moonlet_radius) continue;
			pt.vx 		= 0;
			pt.vy 		= -1.5*pt.x*OMEGA;
			pt.vz 		= 0;
			pt.ax 		= 0; pt.ay 		= 0; pt.az 		= 0;
			double radius 	= tools_powerlaw(particle_radius_min,particle_radius_max,particle_radius_slope);
			pt.r 		= radius;						// m
			double	particle_mass = particle_density*4./3.*M_PI*radius*radius*radius;
			pt.m 		= particle_mass; 	// kg
			particles_add(pt);
			mass += particle_mass;
		}
	}
}

double coefficient_of_restitution_bridges(double v){
	// v is in units of [m/s]
	double eps = 0.32*pow(fabs(v)*100.,-0.234);
	if (eps>1) eps=1;
	if (eps<0) eps=0;
	return eps;
}

double shift_x = 0;
double shift_y = 0;

void output_moonlet(char* filename){
	struct particle moon = particles[0]; // Moonlet
	FILE* f = fopen(filename, "a");
	fprintf(f,"%e\t",t);					// Time
	fprintf(f,"%e\t%e\t%e\t",moon.x+shift_x, moon.y+shift_y, moon.z );	// Position
	fprintf(f,"%e\t%e\t%e\t",moon.vx,moon.vy-1.5*OMEGA*shift_x,moon.vz);	// Velocity
	fprintf(f,"%e\t%e\t%e\t",moon.ax,moon.ay,moon.az);	// Acceleration
	fprintf(f,"\n");
	fclose(f);
}

void problem_inloop(){
}

int output_num_moonlet = 0;
void problem_output(){
	// Recenter moonlet
	struct particle moonlet = particles[0];
	double aO = 2.*moonlet.vy/OMEGA + 4.*moonlet.x;	// Center of epicyclic motion
	double bO = moonlet.y - 2.*moonlet.vx/OMEGA;
	shift_x += aO;
	shift_y += bO;
	for (int i=0;i<N;i++){
		particles[i].x  -= aO;
		particles[i].y  -= bO;
		particles[i].vy += 1.5*OMEGA*aO;
	}

	// Output stuff
	if (output_check(1e-2*2.*M_PI/OMEGA)){
#ifdef MPI
		if (mpi_id==0){
#endif 	// MPI
			output_moonlet("moonlet.txt");
#ifdef MPI
		}
#endif 	// MPI
	}
	if (output_check(1e-2*2.*M_PI/OMEGA)){
		output_timing();
	}
	if (output_check(2.*M_PI/OMEGA)){
		char filename[1024];
		sprintf(filename, "positions_t%05d.txt", output_num_moonlet);
		output_num_moonlet++;
		output_ascii(filename);
	}
}

void problem_finish(){
}

/**
 * @file 	problem.c
 * @brief 	Example problem: self-gravity disc.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	A self-gravitating disc is integrated using
 * the leap frog integrator. This example is also compatible with 
 * the Wisdom Holman integrator. Collisions are not resolved.
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"
#include "communication_mpi.h"
#include "tree.h"
#include "input.h"
#include "tools.h"

double ss_pos[6][3] = 
{
	{-4.06428567034226e-3,	-6.08813756435987e-3,	-1.66162304225834e-6	}, // Sun
	{+3.40546614227466e+0,	+3.62978190075864e+0,	+3.42386261766577e-2	}, // Jupiter
	{+6.60801554403466e+0,	+6.38084674585064e+0,	-1.36145963724542e-1	}, // Saturn
	{+1.11636331405597e+1,	+1.60373479057256e+1,	+3.61783279369958e-1	}, // Uranus
	{-3.01777243405203e+1,	+1.91155314998064e+0,	-1.53887595621042e-1	}, // Neptune
	{-2.13858977531573e+1,	+3.20719104739886e+1,	+2.49245689556096e+0	}  // Pluto
};
double ss_vel[6][3] = 
{
	{+6.69048890636161e-6,	-6.33922479583593e-6,	-3.13202145590767e-9	}, // Sun
	{-5.59797969310664e-3,	+5.51815399480116e-3,	-2.66711392865591e-6	}, // Jupiter
	{-4.17354020307064e-3,	+3.99723751748116e-3,	+1.67206320571441e-5	}, // Saturn
	{-3.25884806151064e-3,	+2.06438412905916e-3,	-2.17699042180559e-5	}, // Uranus
	{-2.17471785045538e-4,	-3.11361111025884e-3,	+3.58344705491441e-5	}, // Neptune
	{-1.76936577252484e-3,	-2.06720938381724e-3,	+6.58091931493844e-4	}  // Pluto
};

double ss_mass[6] =
{
	1.00000597682, 	// Sun + inner planets
	1./1047.355,	// Jupiter
	1./3501.6,	// Saturn
	1./22869.,	// Uranus
	1./19314.,	// Neptune
	0.		// Pluto
};

const double k	 	= 0.01720209895;	// Gaussian constant 
double energy();
double energy_init; 
#ifdef INTEGRATOR_RADAU15
extern double integrator_epsilon;
extern int integrator_force_is_velocitydependent;
#endif

void move_to_heliocentric_frame(){
	// Move to heliocentric frame (required by WHM)
	for (int i=1;i<N;i++){
		particles[i].x -= particles[0].x;	particles[i].y -= particles[0].y;	particles[i].z -= particles[0].z;
		particles[i].vx -= particles[0].vx;	particles[i].vy -= particles[0].vy;	particles[i].vz -= particles[0].vz;
	}
	particles[0].x = 0;	particles[0].y = 0;	particles[0].z = 0;
	particles[0].vx= 0;	particles[0].vy= 0;	particles[0].vz= 0;
}

void problem_init(int argc, char* argv[]){
	// Setup constants
	G		= k*k;

	// Setup homog. sphere
	tmax		= 1e5*365.;
	dt		= 365./input_get_double(argc,argv,"timesteps",100);
#ifdef INTEGRATOR_RADAU15
	integrator_epsilon = input_get_double(argc,argv,"epsilon",0);
	integrator_force_is_velocitydependent = 0;
#endif
	boxsize = 500;	
	init_box();

	// Initial conditions
	for (int i=0;i<6;i++){
		struct particle p;
		p.x  = ss_pos[i][0]; 		p.y  = ss_pos[i][1];	 	p.z  = ss_pos[i][2];
		p.vx = ss_vel[i][0]; 		p.vy = ss_vel[i][1];	 	p.vz = ss_vel[i][2];
		p.ax = 0; 			p.ay = 0; 			p.az = 0;
		p.m  = ss_mass[i];
		particles_add(p); 
	}
	tools_move_to_center_of_momentum();
	energy_init = energy();
#ifdef INTEGRATOR_WH
	move_to_heliocentric_frame();
#endif // INTEGRATOR_WH


	mpf_set_default_prec(512);
}

double energy(){
	mpf_t energy_kinetic;
	mpf_init(energy_kinetic);
	mpf_t energy_potential;
	mpf_init(energy_potential);
	for (int i=0;i<N;i++){
		for (int j=0;j<i;j++){
			double dx = particles[i].x - particles[j].x;
			double dy = particles[i].y - particles[j].y;
			double dz = particles[i].z - particles[j].z;
			double r = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
			mpf_t temp;
			mpf_init(temp);
			mpf_set_d(temp,-G*particles[i].m*particles[j].m/r);
			mpf_add(energy_potential, energy_potential,temp);
			
		}
		double dvx = particles[i].vx;
		double dvy = particles[i].vy;
		double dvz = particles[i].vz;
		mpf_t temp;
		mpf_init(temp);
		mpf_set_d(temp,1./2.*particles[i].m * (dvx*dvx + dvy*dvy + dvz*dvz));
		mpf_add(energy_kinetic, energy_kinetic,temp);
	}
	mpf_add(energy_kinetic,energy_kinetic,energy_potential);

	return mpf_get_d(energy_kinetic);
}
void problem_inloop(){
}

FILE* ofe=NULL;
//FILE* ofe = NULL; 
void problem_finish(){
	printf("\nFinished. Time: %e. Difference to exact time: %e. Timestep: %e\n",t,t-tmax,dt);
#ifdef INTEGRATOR_LEAPFROG
	FILE* of = fopen("energy_leapfrog.txt","a+"); 
#endif
#ifdef INTEGRATOR_RADAU15
	FILE* of = fopen("energy_radau15.txt","a+"); 
#endif
#ifdef INTEGRATOR_WH
	FILE* of = fopen("energy_wh.txt","a+"); 
#endif
	fprintf(of,"%e\t",dt);
	tools_move_to_center_of_momentum();
	fprintf(of,"%e\t",fabs((energy()-energy_init)/energy_init));
		struct timeval tim;
		gettimeofday(&tim, NULL);
		double timing_final = tim.tv_sec+(tim.tv_usec/1000000.0);
	fprintf(of,"%e\t",timing_final-timing_initial);
#ifdef INTEGRATOR_RADAU15
	fprintf(of,"%e\t",integrator_epsilon);
#endif
	fprintf(of,"\n");
	fclose(of);
//	fclose(ofe);
}
void problem_output(){
//	if(output_check(tmax/10000.)){
//		output_timing();
//		if (ofe==NULL){
//			ofe= fopen("energy.txt","a+"); 
//		}
//		fprintf(ofe,"%e\t",t);
//		fprintf(ofe,"%e\t",dt);
//		fprintf(ofe,"%e\t",integrator_epsilon);
//		fprintf(ofe,"%e\t",fabs((energy()-energy_init)/energy_init));
//		fprintf(ofe,"\n");
//	}
}


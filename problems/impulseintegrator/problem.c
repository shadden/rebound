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

extern int Nmax;
int writeBest; 
double calculate_energy(struct particle* _particles, int _N);
void output_statistics();
extern double opening_angle2; /**< Square of the cell opening angle \f$ \theta \f$. */

void problem_init(int argc, char* argv[]){
	// Setup constants
	G 		= 1;		
	dt		= input_get_double(argc,argv,"dt",1./48.);
	writeBest	= input_get_int(argc,argv,"write",0);

	
	int _N 		= input_get_int(argc,argv,"N",4096);
	double M 	= 1;
	double R 	= 3./64.*M_PI;  // --> Energy = 1
	tmax		= 8;
	boxsize 	= 200.*R;
	softening 	= 0.025;		
#ifdef TREE
	opening_angle2	= 1;
#endif

	init_box();
	
	
	
	tools_init_plummer(_N, M, R);
	tools_move_to_center_of_momentum();

	system("rm -v ascii.txt");
	system("rm -v E4E8.txt");
}

void problem_inloop(){
}
double calculate_mininteractiontime(){
	double mininteractiontime = 1e12;
	for (int i=0;i<N;i++){
		for (int j=0;j<i;j++){
			double dx = particles[i].x - particles[j].x;
			double dy = particles[i].y - particles[j].y;
			double dz = particles[i].z - particles[j].z;
			double r = sqrt(dx*dx + dy*dy + dz*dz);
			double dvx = particles[i].vx - particles[j].vx;
			double dvy = particles[i].vy - particles[j].vy;
			double dvz = particles[i].vz - particles[j].vz;
			double v = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
			double interactiontime = fabs(r/v);
			if (interactiontime<mininteractiontime){
				mininteractiontime = interactiontime;
			}
		}
	}
	return mininteractiontime;
}

struct vec3 {
	double x;
	double y;
	double z;
};

void output_statistics(){
	printf("\n\n");
	printf("\t Number of particles: \t%d\n",N);
	double mass = 0;
	for (int i=0;i<N;i++){
		mass += particles[i].m;
	}
	printf("\t Total mass:           \t%e\n",mass);
	// Algorithm with reduced roundoff errors (see wikipedia)
	struct vec3 A = {.x=0, .y=0, .z=0};
	struct vec3 Q = {.x=0, .y=0, .z=0};
	for (int i=0;i<N;i++){
		struct vec3 Aim1 = A;
		struct particle p = particles[i];
		A.x = A.x + (p.vx-A.x)/(double)(i+1);
		A.y = A.y + (p.vy-A.y)/(double)(i+1);
		A.z = A.z + (p.vz-A.z)/(double)(i+1);
		Q.x = Q.x + (p.vx-Aim1.x)*(p.vx-A.x);
		Q.y = Q.y + (p.vy-Aim1.y)*(p.vy-A.y);
		Q.z = Q.z + (p.vz-Aim1.z)*(p.vz-A.z);
	}
	Q.x = sqrt(Q.x/(double)N);
	Q.y = sqrt(Q.y/(double)N);
	Q.z = sqrt(Q.z/(double)N);
	double totQ = sqrt(Q.x*Q.x + Q.y*Q.y + Q.z*Q.z);
	printf("\t Velocity dispersion: \t%e\n",totQ);

	double energy_kinetic = 0;
	double energy_potential = 0;
#pragma omp parallel for reduction(+:energy_potential) reduction(+:energy_kinetic) schedule(guided)
	for (int i=0;i<N;i++){
		for (int j=0;j<i;j++){
			double dx = particles[i].x - particles[j].x;
			double dy = particles[i].y - particles[j].y;
			double dz = particles[i].z - particles[j].z;
			double r = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
			energy_potential += -G*particles[i].m*particles[j].m/r;
		}
		double dvx = particles[i].vx;
		double dvy = particles[i].vy;
		double dvz = particles[i].vz;
		energy_kinetic += 1./2. * particles[i].m* (dvx*dvx + dvy*dvy + dvz*dvz);
	}
	
	printf("\t Energy (kinetic): \t%e\n",energy_kinetic);
	printf("\t Energy (potential): \t%e\n",energy_potential);
	printf("\t Energy (total): \t%e\n",energy_potential+energy_kinetic);
	printf("\n\n");
	
	FILE* of = fopen("ascii.txt","a+"); 
	fprintf(of,"%e\t",t);
	fprintf(of,"%e\t",energy_kinetic);
	fprintf(of,"%e\t",energy_potential);
	fprintf(of,"\n");
	fclose(of);
	
}

int compare (const void * a, const void * b){
	if (*(double*)a < *(double*)b ) return -1;
	if (*(double*)a > *(double*)b ) return 1;
	return 0;
}

void output_radii(){
	tools_move_to_center_of_momentum();
	double* radii = malloc(sizeof(double)*N);
	for (int i=0;i<N;i++){
		double r = sqrt(particles[i].x*particles[i].x + particles[i].y*particles[i].y + particles[i].z*particles[i].z);
		radii[i] = r;
	}
	qsort (radii, N, sizeof(double), compare);
	
	FILE* of = fopen("radii.txt","a+"); 
	fprintf(of,"%e\t",t);
	int index10 = (int)ceil(0.1*(double)N);
	if (index10>=N) index10 = N-1;
	int index50 = (int)ceil(0.5*(double)N);
	if (index50>=N) index50 = N-1;
	int index90 = (int)ceil(0.9*(double)N);
	if (index90>=N) index90 = N-1;
	fprintf(of,"%e\t",radii[index10]);
	fprintf(of,"%e\t",radii[index50]);
	fprintf(of,"%e\t",radii[index90]);
	
//	for (int i=0;i<N;i++){
//		fprintf(of,"%e\t%d\t%e\n",t/t0,i,radii[i]/r0);
//	}

	fprintf(of,"\n");
	fclose(of);
}


struct vec3o {
	double x;
	double y;
	double z;
};

void output_E4E8(){
	FILE* of = fopen("E4E8.txt","a+"); 
	for (int i=0;i<N;i++){
		double energy_potential = 0;
		double energy_kinetic = 0;
		for (int j=0;j<N;j++){
			if (i==j) continue;
			double dx = particles[i].x - particles[j].x;
			double dy = particles[i].y - particles[j].y;
			double dz = particles[i].z - particles[j].z;
			double r = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
			energy_potential += -G*particles[i].m*particles[j].m/r;
		}
		double dvx = particles[i].vx;
		double dvy = particles[i].vy;
		double dvz = particles[i].vz;
		energy_kinetic += 1./2. * particles[i].m* (dvx*dvx + dvy*dvy + dvz*dvz);
		fprintf(of,"%e\t%e\n",particles[i].E4,energy_kinetic + energy_potential);
	}
	fclose(of);
}

void problem_finish(){
	printf("\nFinished. Time: %e. Difference to exact time: %e. Timestep: %e\n",t,t-tmax,dt);
	if (writeBest){
		output_binary("best.bin");
	}
	output_E4E8();
}

void set_E4(){
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		double energy_potential = 0;
		double energy_kinetic = 0;
		for (int j=0;j<N;j++){
			if (i==j) continue;
			double dx = particles[i].x - particles[j].x;
			double dy = particles[i].y - particles[j].y;
			double dz = particles[i].z - particles[j].z;
			double r = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
			energy_potential += -G*particles[i].m*particles[j].m/r;
		}
		double dvx = particles[i].vx;
		double dvy = particles[i].vy;
		double dvz = particles[i].vz;
		energy_kinetic += 1./2. * particles[i].m* (dvx*dvx + dvy*dvy + dvz*dvz);
		particles[i].E4 = energy_kinetic + energy_potential;
	}
}

void problem_output(){
	//printf("Min interaction time: %e\n",calculate_mininteractiontime());
	output_timing();
	if (t<=4.&& t+dt>4.){
		set_E4();
	}
	if (output_check(tmax)) output_statistics();
//#	if (output_check(tmax/1000.)) output_radii();
}


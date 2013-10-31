/**
 * @file 	problem.c
 * @brief 	Example problem: circular orbit.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example uses the Wisdom Holman integrator
 * to integrate particles on a circular orbit in a fixed 
 * potential.
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
#include <string.h>
#include "main.h"
#include "problem.h"
#include "input.h"
#include "tools.h"
#include "output.h"
#include "particle.h"
#include "boundaries.h"
#include "integrator.h"

// Star and planet (note, those wont be updated after they have been inserted)
void additional_forces();
const double Msaturn	= 0.00028588598; 	// Mass of Saturn in solar masses 
const double Rsaturn	= 0.00038925688; 	// Radius of Saturn in AU
const double Rsun   	= 0.0046524726; 	// Radius of Sun in AU
const double J2saturn	= 16298e-6; 		// J2 of Saturn (Murray and Dermott p 531) 
const double Rring_inner_saturn	= 0.00049905791;// Inner Saturn ring edge (C) in AU from wikipedia
const double Rring_outer_saturn	= 0.00091431783;// Outer Saturn ring edge (outside Keeler gap) in AU from wikipedia
const double Rhill_saturn	= 0.43767071; 	// Hill sphere of Saturn in AU

double Mplanet;
double aplanet;
double Pplanet;
double Rplanet;
double Rhill_planet;
double Rring_inner_planet;
double Rring_outer_planet;
double Obliquity;
double Rstar;
double J2planet;
double betaparticles;

extern double integrator_epsilon;

void problem_init(int argc, char* argv[]){
	// Setup constants
	aplanet				= input_get_double(argc,argv,"aplanet",0.2); 
	integrator_epsilon		= input_get_double(argc,argv,"integrator_epsilon",1e-2); 
	boxsize 			= 5.*aplanet;
	N_active      			= 2;
	problem_additional_forces 	= additional_forces;
	init_box();

	// Initial conditions
	// Star
	Rstar = input_get_double(argc,argv,"Rstar",1)*Rsun;
	struct particle star;
	star.x  = 0; star.y  = 0; star.z  = 0;
	star.vx = 0; star.vy = 0; star.vz = 0;
	star.ax = 0; star.ay = 0; star.az = 0;
	star.m  = input_get_double(argc,argv,"Mstar",1);
	particles_add(star);


	// Planet 
	struct particle planet;
	Mplanet = input_get_double(argc,argv,"Mplanet",1)*Msaturn;
	planet.m  = Mplanet;
	planet.x  = aplanet;
	planet.y  = 0.; 
	planet.z  = 0.;
	planet.ax = 0; planet.ay = 0; planet.az = 0;
	planet.vx = 0;
	planet.vy = sqrt(G*(star.m+planet.m)/planet.x);
	planet.vz = 0;
	particles_add(planet);
	Rhill_planet 		= aplanet*pow(planet.m/star.m/3.,1./3.);
	Pplanet			= 2.*M_PI*sqrt(aplanet*aplanet*aplanet/(1.+Mplanet)/G);
	Rplanet 		= input_get_double(argc,argv,"Rplanet",1)*Rsaturn;
	J2planet 		= input_get_double(argc,argv,"J2planet",1)*J2saturn;
	Obliquity 		= input_get_double(argc,argv,"Obliquity",0)/180.*M_PI;

	// Ring particles
	betaparticles 		= input_get_double(argc,argv,"betaparticles",0);
	Rring_inner_planet 	= input_get_double(argc,argv,"Rring_inner_planet",1)*Rring_inner_saturn;
	Rring_outer_planet 	= input_get_double(argc,argv,"Rring_outer_planet",1)*Rring_outer_saturn;
	int _N = input_get_int(argc,argv,"N",10);

	for(int i=0;i<_N;i++){
		struct particle tp = planet;
		tp.m  = 0;
		double r 	= Rring_inner_planet + (Rring_outer_planet-Rring_inner_planet)*i/(_N-1); 
		if (i!=0&&i!=_N-1){
			r += 0.1*tools_normal(1)*(Rring_outer_planet-Rring_inner_planet)/(_N-1);
		}
		double v	= sqrt(G*planet.m / r *(1.+3.*J2planet*Rplanet*Rplanet/2./r)); // correction puts particles on circular orbits when J2 is present
		double phi 	= tools_uniform(0,2.*M_PI);
		// Circular orbit around planet, randomize azimuth
		double x  =  r * cos(phi);
		double y  =  r * sin(phi);
		double vx = -v * sin(phi);
		double vy =  v * cos(phi);
		// Rotate around y axis by Obliquity
		tp.x	+=  x  * cos(Obliquity);
		tp.y  	+=  y;
		tp.z	+= -x  * sin(Obliquity);
		tp.vx	+=  vx * cos(Obliquity);
		tp.vy  	+=  vy;
		tp.vz	+= -vx * sin(Obliquity);
		particles_add(tp);
	}

	dt = 2.*M_PI*sqrt(Rring_inner_planet*Rring_inner_planet*Rring_inner_planet/Mplanet/G)*input_get_double(argc,argv,"dtfrac",1e-1);
	tmax = 1e5*2.0*M_PI*sqrt(aplanet*aplanet*aplanet/G/star.m);
	tools_move_to_center_of_momentum();

	output_prepare_directory();
}
		
void force_J2(){
	if (J2planet==0) return;
	// Star 
	const struct particle planet = particles[1];				// cache
#pragma omp parallel for
	for (int i=2;i<N;i++){
		const struct particle p = particles[i]; 			// cache
		const double sprx  = p.x-planet.x;
		const double spry  = p.y-planet.y;
		const double sprz  = p.z-planet.z;
		const double prx  = sprx*cos(-Obliquity) + sprz*sin(-Obliquity);
		const double pry  = spry;
		const double prz  =-sprx*sin(-Obliquity) + sprz*cos(-Obliquity);
		const double pr2   = prx*prx + pry*pry + prz*prz; 		// distance^2 relative to planet
		const double fac  = 3.*G*J2planet*planet.m*Rplanet*Rplanet/2./pow(pr2,3.5);

		const double pax = fac*prx*(prx*prx + pry*pry - 4.*prz*prz);
		const double pay = fac*pry*(prx*prx + pry*pry - 4.*prz*prz);
		const double paz = fac*prz*(3.*(prx*prx + pry*pry) - 2.*prz*prz);
		
		particles[i].ax += pax*cos(Obliquity) + paz*sin(Obliquity);
		particles[i].ay += pay;
		particles[i].az +=-pax*sin(Obliquity) + paz*cos(Obliquity);
	}
}


void force_radiation(){
	if (betaparticles==0) return;
	// Star 
	const struct particle star = particles[0];				// cache
//#define SHADOW
#ifdef SHADOW
	const struct particle planet = particles[1];				// cache
#endif // SHADOW
#pragma omp parallel for
	for (int i=2;i<N;i++){
		const struct particle p = particles[i]; 			// cache
		const double prx  = p.x-star.x;
		const double pry  = p.y-star.y;
		const double prz  = p.z-star.z;
		const double pr   = sqrt(prx*prx + pry*pry + prz*prz); 	// distance relative to star

#ifdef SHADOW
		const double plrx = planet.x-star.x;
		const double plry = planet.y-star.y;
		const double plrz = planet.z-star.z;
		const double plr  = sqrt(plrx*plrx + plry*plry + plrz*plrz); 	// distance of planet relative to star
		int IN_SHADOW = 0;
		// Find out if particle is in the shadow of the planet
		// i.e., find out if the radial vector to the particle passes through the planet
		// If the particle isn't as far as the planet's center, it's definitely not in shadow.
		if (pr < plr){
			IN_SHADOW = 0;
		}else{
			// If the angle between the particle's position vector relative to the star and the center of
			// the planet's position vector relative to the center of the star is bigger than Rplanet/pr, then
			// the planet is not in shadow.  Otherwise, it is.
			const double max_angle_from_planet_center = Rplanet / pr; // planet radius / distance from star
			// Call vector to planet center "Pl" and vector to particle "p"
			// angle from planet center is arccos[Pl dot p / (plr * pr)]
			const double dotprod = prx*plrx + pry*plry + prz*plrz;
			const double cos_angle_from_planet_center = dotprod / (plr*pr);
			const double angle_from_planet_center = acos(cos_angle_from_planet_center);
			if (angle_from_planet_center < max_angle_from_planet_center)
				IN_SHADOW = 1;
		}
		if (IN_SHADOW == 1) continue;
#endif // SHADOW

		// No shadow; calculate force modification:
		const double prvx = p.vx-star.vx;
		const double prvy = p.vy-star.vy;
		const double prvz = p.vz-star.vz;

		//const double beta	= 0.114;    				// appropriate for 5-micron dust (Burns et al. 1979, Eq 19, Q=rho=1)
		const double c 		= 10064.915; 				// speed of light.
		const double rdot 	= (prvx*prx + prvy*pry + prvz*prz)/pr; 	// radial velocity relative to star
		const double F_r 	= betaparticles*G*star.m/(pr*pr);

		// Equation (5) of Burns, Lamy, Soter (1979)
		particles[i].ax += F_r*((1.-rdot/c)*prx/pr - prvx/c);
		particles[i].ay += F_r*((1.-rdot/c)*pry/pr - prvy/c);
		particles[i].az += F_r*((1.-rdot/c)*prz/pr - prvz/c);

	}
}

void additional_forces(){
	force_J2();
	force_radiation();
}
						
void output_append_ascii_rings(char* filename){
	FILE* of = fopen(filename,"a"); 
	struct particle planet = particles[1];
	for (int i=2;i<N;i++){
		struct particle p = particles[i];
		
		fprintf(of,"%e\t%e\t%e\t%e\n",t,p.x-planet.x,p.y-planet.y,p.z-planet.z);
	}
	fclose(of);
}

int burst_count = -1;
int burst = 0;
int burst_max = 1000;
void problem_output(){
	if(output_check(2.*M_PI)){
		const struct particle star 	= particles[0];				// cache
		const struct particle planet 	= particles[1];				// cache
		for(int i=2;i<N;i++){
			const struct particle p = particles[i]; 			// cache
			const double psrx  = p.x-star.x;
			const double psry  = p.y-star.y;
			const double psrz  = p.z-star.z;
			const double pprx  = p.x-planet.x;
			const double ppry  = p.y-planet.y;
			const double pprz  = p.z-planet.z;
			const double psr  = sqrt(psrx*psrx + psry*psry + psrz*psrz); 	// distance of particle relative to star
			const double ppr  = sqrt(pprx*pprx + ppry*ppry + pprz*pprz); 	// distance of particle relative to planet
			if (ppr<Rplanet || psr<Rstar){
				// Remove particle
				particles[i] = particles[N-1];
				i--;
				N--;
			}
		}
	}
	
	if(N<3){
		printf("Last testparticle lost.\nExiting.\n");
		exit(0);
	}
	// Output stuff
	if(output_check(2.*M_PI)){
		output_timing();
	}
	if(output_check(100.*Pplanet)){
		burst = burst_max;
		burst_count++;
	}
	if (burst){
		char filename[128];
		sprintf(filename,"ascii_%04d.txt",burst_count);
		output_append_ascii_rings(filename);
		burst--;
	}
	if (fabs(particles[0].x)>boxsize_x/10. || fabs(particles[0].y)>boxsize_y/10. || fabs(particles[0].z)>boxsize_z/10. ){
		tools_move_to_center_of_momentum();
	}
}
void problem_inloop(){
}

void problem_finish(){
}

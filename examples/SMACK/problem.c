/**
 * @file 	problem.c
 * @brief 	Example problem: shearing sheet.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 		Akihiko Fujii <akihiko.fujii@nao.ac.jp>
 * @detail 	This problem uses shearing sheet boundary
 * conditions. Particle properties resemble those found in 
 * Saturn's rings. The collision resolve method has been changed
 * such that particle clumps behave more realistically.
 * 
 * @section 	LICENSE
 * Copyright (c) 2014 Hanno Rein, Shangfei Liu, Akihiko Fujii
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
#include "collision_resolve.h"
#include "output.h"
#include "communication_mpi.h"
#include "tree.h"
#include "tools.h"
#include "display.h"

extern int Nmax;

#define numbins 21
#define extra 30

const int nowave = 1; // Flag to remove size distribution wave by extrapolation
const double logbinsize = 0.1; // bin size for size distribution
const int shift = 3;		// for building size distribution
const int rho = 3.0E3;		// Density (kg m^-3)

double sizebins[numbins];		// size bins
double bigbins[numbins+extra];	// extrapolated size bins
double mass[numbins];			// mass bins
double bigmass[numbins+extra];	// extrapolated mass bins
double Qd[numbins+extra];		// minimum projectile kinetic energy
double Ecol[numbins+extra][numbins+extra];  // collision energy
double Qsuper[numbins+extra];// threshold for supercatastrophic fragmentation
double Mlr[numbins+extra][numbins+extra]; // mass of largest fragment
double powbigbins[numbins+extra]; // used to normalize fragment distribution
double powsizebins[numbins];	// used to normalize fragment distribution
double fragdistsum[numbins];	// used to normalize fragment distribution

#include "smack.c"

void problem_init(int argc, char* argv[]){
	// Setup constants
	N_active = 1;		// Number of active particles
	N_tree_fixed = 1;	// Number of fixed particles in tree
	boxsize = 180;		// AU
	root_nx = 2;
	root_ny = 2;
	dt = 1.0/0.159;		// Timestep (years->code units)
	tmax = (1.0e2+1)*dt;// Max time
	init_box();
	
	// Change collision_resolve routing from default.
	collision_resolve = collision_resolve_single_fragment;
	
	// Initial conditions for star
	struct particle star;
	star.x 		= 0; star.y 	= 0; star.z	= 0;
	star.vx 	= 0; star.vy 	= 0; star.vz 	= 0;
	star.ax 	= 0; star.ay 	= 0; star.az 	= 0;
	star.m 		= 1.0;
	star.r		= 1.0*0.0046491;	// Rsun

	// Add the star to the particle array, but not to the tree
	Nmax += 128;
	particles = realloc(particles,sizeof(struct particle)*Nmax);
	particles[N] = star;
	N++;
	
	// Orbital parameters of disk
	double amin = 50.0;			// Semi-major axis
	double amax = 60.0;
	double emin = 0.0;			// Eccentricity
	double emax = 0.1;	
	double imin = 0.0;			// Inclination
	double imax = emax/2.; 
	double Orange = 2*M_PI;		// Longitude of ascending node
	double orange = 2*M_PI;		// Argument of periapse
	double Mrange = 2*M_PI;		// Mean anomaly
	double a,e,i,OMEGA,omega,M;
	
	// Superparticle parameters
	int numswarms = 1000;		// Number of superparticles
	double swarmr = 1.0;		// Superparticle radius
	double height = 2.77137;	// Height of disk
	double factor = 3*height/(4*swarmr);	// Filling factor
	double initOD = 1.0e-2;					// Initial optical depth
	double Aeach = initOD*M_PI*swarmr*swarmr*1.5e11*1.5e11/factor; // Cross-sectional area 
	double p0 = 2.5; // Initial size distribution index
	double distsum;
	double C;
	double sdist[numbins];
	double smass;	
	
	// Create size bin arrays
	for (int i=0;i<numbins;i++){
		sizebins[i] = pow(10,i*logbinsize-shift);
	}
	if (nowave) { // Extrapolate to smaller sizes to remove size distribution wave
		for (int i=0; i<numbins+extra; i++) {
			bigbins[i] = pow(10,i*logbinsize-shift-extra*logbinsize);
		}
	}
	
	// Collision energy pre-calculations
	const double Scon = 3.5e3;	// Strength regime coefficient
	const double Gcon = 3.0e-8;	// Gravity regime coefficient
	const double spow = -0.38;	// Strength regime exponent
	const double gpow = 1.36;	// Gravity regime exponent
	const double fKE = 1.0;		// Fraction of kinetic energy used
	{
		double bigvolume[numbins+extra];	// extrapolated volume bins
		for (int i=0; i<numbins+extra; i++) {
			if (i>=extra) {
				mass[i-extra] = rho*(4./3)*M_PI*pow(sizebins[i-extra]/2.,3); // small array of mass bins
			}
			bigvolume[i] = (4./3)*M_PI*pow(bigbins[i]/2.,3); // large array of volume bins
			bigmass[i] = rho*bigvolume[i]; // large array of mass bins
			double Qs = Scon*pow(bigbins[i]*100./2,spow); // calculate internal binding energy
			double Qg = Gcon*rho*pow(bigbins[i]*100./2,gpow); // calculate gravitational binding energy
			Qd[i] = (1./fKE)*(Qs+Qg); // calculate minimum projectile kinetic energy
			Qsuper[i] = -2.0*Qd[i]*(0.1-1);
		}
	}
	
	// Fragment distribution pre-caculations
	const double n0 = 2.8;			// Index of fragment distribution
	for (int j = 0; j < numbins+extra; j++) {
		powbigbins[j] = pow(bigbins[j],3.-n0);
	}
	for (int j = 0; j < numbins; j++) {
		powsizebins[j] = pow(sizebins[j],-1.*n0);
	}
	
	for (int lr = 0;  lr < numbins; lr++) {
		for (int k = 0; k <= (lr+extra); k++) {
			fragdistsum[lr] += powbigbins[k];
		}	
	}
	
	// Setup particle structures
#ifdef MPI
	for(int N=0;N<numswarms/mpi_num;N++) {
#else
	for(int N=0;N<numswarms;N++) {
#endif
		// Create orbital element distributions
		a = tools_uniform(amin,amax);
		e = tools_uniform(emin,emax);
		i = tools_uniform(imin,imax);
		OMEGA = tools_uniform(0,Orange);
		omega = tools_uniform(0,orange);
		M = tools_uniform(0,Mrange);
		
		// Fill size distributions
		smass = 0.0;	// Superparticle mass
		distsum = 0.0;
		for (int j=0; j<numbins; j++) {
			distsum += pow(bigbins[j],-1*p0)*pow(bigbins[j],2);
		}
		C = (4*Aeach)/(M_PI*distsum);	// Normalize size distribution
		for (int i=0; i<numbins; i++) {
			sdist[i] = C*pow(sizebins[i],-p0);
			smass += (M_PI/6)*rho*pow(sizebins[i],3)*sdist[i];
		}
		
		smass /= 1.98892e30; // Convert superparticles mass to solar units
		struct particle p;
		p = tools_orbit2p(a,e,i,omega,OMEGA,M,star.m,smass); // Calculate Cartesian coords
		p.m = smass;	// Superparticle mass
		p.r = swarmr;	// Superparticle radius
		p.number = N + N_active;	// Superparticle number
#ifdef MPI
		p.number = N + N_active + mpi_id*numswarms/mpi_num;
#endif
		p.ncol = 0;	// Number of collisions
		for (int i=0; i<numbins; i++) {
			p.sdist[i] = sdist[i];	// Superparticle size distribution
		}	
		particles_add(p);
		
#ifdef MPI			
		communication_mpi_distribute_particles();
#endif
	}
	
}


void problem_inloop(){
}

void problem_output(){
	output_timing();
	if (output_check(dt*100.0)) {
		output_smack_sizedist("sizedist.txt");
	}
}

void problem_finish(){
}

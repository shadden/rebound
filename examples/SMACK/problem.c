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

extern long collisions_Nlog;
struct line {
	double slope;
	double intercept;
};

/**
 * This function fits a set of numbers to a straight line y=ax+b, returning a and b
 * @author: Erika Nesvold
 * @param x[] x data
 * @param y[] y data
 * @param size Length of data
 * @return line
 */
struct line tools_linefit(double x[], double y[], int size);

/**
 * This function calculates x,y,z,vx,vy,vz given a,e,i,omega,OMEGA,M
 */

struct particle tools_orbit2p(float a, float e, float i, float omega, float OMEGA, float M, float Ms, float Mp);

/**
 * Appends the size distributions of the swarms to an ASCII file.
 * @author: Erika Nesvold
 * @param filename Output filename
 */
void output_sizedist_append(char* filename);

void collision_resolve_single_fragment(struct collision c);

extern double opening_angle2;

int numbins = 21;
float logbinsize = 0.1;
int shift = 3;

void problem_init(int argc, char* argv[]){
	// Setup constants
	N_active = 1;
#ifndef COLLISIONS_NONE	
	N_tree_fixed = 1;
#endif
	boxsize = 180;	// AU
	root_nx = 2;
	root_ny = 2;
	dt = 1.0/0.159; // Timestep (years->code units)
	tmax = (1.0e2+1)*dt;
	init_box();
	
	// Change collision_resolve routing from default.
	collision_resolve = collision_resolve_single_fragment;
	
	// Initial conditions for star
	struct particle star;
	star.x 		= 0; star.y 	= 0; star.z	= 0;
	star.vx 	= 0; star.vy 	= 0; star.vz 	= 0;
	star.ax 	= 0; star.ay 	= 0; star.az 	= 0;
	star.m 		= 1.0;
#ifndef COLLISIONS_NONE
	star.r		= 1.0*0.0046491;	// Rsun
#endif
	particles_add_local_not_tree(star);
	
	// Collision parameters
	float p0 = 2.5;
	float rho = 3000.;
	int nowave = 1;
	int extra = 30;
	
	// Orbital parameters
	float amin = 50.0;
	float amax = 60.0;
	float emin = 0.0;
	float emax = 0.1;
	float imin = 0.0;
	float imax = emax/2.; 
	float Orange = 2*M_PI;	
	float orange = 2*M_PI;
	float Mrange = 2*M_PI;
	double a,e,i,OMEGA,omega,M;

	// Swarm parameters
	int numswarms = 100;
	//float swarmr = pow(10.0,-(1./3));	// AU
	float swarmr = 1.0;
	float height = 2.77137; // ecc 0.2
	float factor = 3*height/(4*swarmr);
	float initOD = 1.0e-2; // Initial optical depth
	float Aeach = initOD*M_PI*swarmr*swarmr*1.5e11*1.5e11/factor;  
	float distsum;
	float C;
	double sdist[numbins];
	double smass;
	
	// Size bin array
	float sizebins[numbins];
	float bigbins[numbins+extra];
	for (int i=0;i<numbins;i++){
		sizebins[i] = pow(10,i*logbinsize-shift);
	}
	if (nowave) {
		for (int i=0; i<numbins+extra; i++) {
			bigbins[i] = pow(10,i*logbinsize-shift-extra*logbinsize);
		}
	}
	
	
	// Setup particle structures
#ifdef MPI
	for(int N=0;N<numswarms/mpi_num;N++) {
#else
		for(int N=0;N<numswarms;N++) {
#endif
	        a = tools_uniform(amin,amax);
			e = tools_uniform(emin,emax);
			i = tools_uniform(imin,imax);
			OMEGA = tools_uniform(0,Orange);
			omega = tools_uniform(0,orange);
			M = tools_uniform(0,Mrange);
			
			// Size distributions
			smass = 0.0;
			distsum = 0.0;
			for (int j=0; j<numbins; j++) {
				distsum += pow(bigbins[j],-1*p0)*pow(bigbins[j],2);
			}
			C = (4*Aeach)/(M_PI*distsum);
			for (int i=0; i<numbins; i++) {
				sdist[i] = C*pow(sizebins[i],-p0);
				smass += (M_PI/6)*rho*pow(sizebins[i],3)*sdist[i];
			}
			
			smass /= 1.98892e30; // convert to solar units
			struct particle p;
			p = tools_orbit2p(a,e,i,omega,OMEGA,M,star.m,smass);
			p.m = smass;
#ifndef COLLISIONS_NONE
			p.r = swarmr;
			p.number = N + N_active;
#ifdef MPI
			p.number = N + N_active + mpi_id*numswarms/mpi_num;
#endif
			p.ncol = 0;
			for (int i=0; i<numbins; i++) {
				p.sdist[i] = sdist[i];
			}
#endif		
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
	if (output_check(dt*10.0)) {
		output_sizedist_append("sizedist.txt");
	}
}

void problem_finish(){
}

void output_sizedist_append(char* filename){
#ifndef COLLISIONS_NONE
#ifdef MPI
		char filename_mpi[1024];
		sprintf(filename_mpi,"%s_%d",filename,mpi_id);
		FILE* of = fopen(filename_mpi,"a");
#else // MPI
		FILE* of = fopen(filename,"a");
#endif //MPI
		for (int i=N_active; i<N; i++) {
			fprintf(of,"%e\t%d\t%d\t",t,(int)particles[i].number,(int)particles[i].ncol);
			for (int j=0; j<numbins; j++) {
		  		fprintf(of,"%e\t",particles[i].sdist[j]);
			}
			fprintf(of,"\n");
		}
		fclose(of);
#endif
}

struct line tools_linefit(double x[], double y[], int size){
	
	struct line l;
	
	double totalx = 0.0;
	double totaly = 0.0;
	double totalxy = 0.0;
	double totalxx = 0.0;
	for (int i=0; i<size; i++) {
		totalx += log10(x[i]);
		totaly += log10(y[i]);
		totalxy += log10(x[i])*log10(y[i]);
		totalxx += log10(x[i])*log10(x[i]);
	}
	double avgx = totalx/size;
	double avgy = totaly/size;
	double avgxy = totalxy/size;
	double avgxx = totalxx/size;
	
	double numerator = avgxy - (avgx*avgy);
	double denominator = avgxx - (avgx*avgx);
	l.slope = numerator/denominator;
	l.intercept = avgy-(l.slope*avgx);
	
	return l;
}

		struct particle tools_orbit2p(float a, float e, float i, float omega, float OMEGA, float M, float Ms, float Mp){
			
			struct particle p;
			// Calculate eccentric anomaly
			float error = 1.0E-6;
			float E;
			if (M<M_PI) {
				E = M + e/2.;
			} else {
				E = M - e/2.;
			}
			float ratio = 1.0;
			while (fabs(ratio)>error) {
				ratio = (E-e*sin(E)-M)/(1-e*cos(E));
				E = E - ratio;
			}
			
			// Calculate x and y before rotation
			float x;
			float y;
			x = a*(cos(E)-e);
			y = a*sqrt(1-pow(e,2))*sin(E);
			
			// Rotate through angles
			float P1[3][3] = {{cos(omega),-1*sin(omega),0.},{sin(omega),cos(omega),0.},{0.,0.,1.}};
			float P2[3][3] = {{1.,0.,0.},{0.,cos(i),-1*sin(i)},{0.,sin(i),cos(i)}};
			float P3[3][3] = {{cos(OMEGA),-1*sin(OMEGA),0.},{sin(OMEGA),cos(OMEGA),0.},{0.,0.,1.}};
			float P4[3][3];
			float Q[3][3];
			for (int i=0; i<3; i++) {
				for (int j=0; j<3; j++) {
					P4[i][j]=P3[i][0]*P2[0][j]+P3[i][1]*P2[1][j]+P3[i][2]*P2[2][j];
				}
			}
			for (int i=0; i<3; i++) {
				for (int j=0; j<3; j++) {
					Q[i][j]=P4[i][0]*P1[0][j]+P4[i][1]*P1[1][j]+P4[i][2]*P1[2][j];
				}
			}
			float R[3][1];
			for (int i=0; i<3; i++) {
				R[i][0] = Q[i][0]*x+Q[i][1]*y+Q[i][2]*0.;
			}
			p.x = R[0][0];
			p.y = R[1][0];
			p.z = R[2][0];
			
			// Calculate velocities
			float mu = Ms+Mp;
			float xdot = -1*sqrt(mu*a)*sin(E)/sqrt(x*x+y*y);
			float ydot = sqrt(mu*a)*sqrt(1-e*e)*cos(E)/sqrt(x*x+y*y);
			float V[3][1];
			for (int i=0; i<3; i++) {
				V[i][0] = Q[i][0]*xdot+Q[i][1]*ydot+Q[i][2]*0.;
			}
			p.vx = V[0][0];
			p.vy = V[1][0];
			p.vz = V[2][0];
			
			p.m = 0;
			p.ax = 0;
			p.ay = 0;
			p.az = 0;
#ifndef COLLISIONS_NONE
			p.r = 0;
			p.lastcollision=0;
#endif
			
			return p;
		}
		

double mfp = 0.0; /* Average mean free path */

void collision_resolve_single_fragment(struct collision c){
#ifndef COLLISIONS_NONE
	if (c.p1*c.p2 == 0) {
		return;	// Don't collide with star
	}
	if (c.p1 == 1 || c.p2 == 1) {
		return;	// Don't collide with planet
	}
	struct particle p1 = particles[c.p1];
	struct particle p2;
#ifdef MPI
	int isloc = communication_mpi_rootbox_is_local(c.ri);
	if (isloc==1){
#endif // MPI
		p2 = particles[c.p2];
#ifdef MPI
	}else{
		int root_n_per_node = root_n/mpi_num;
		int proc_id = c.ri/root_n_per_node;
		p2 = particles_recv[proc_id][c.p2];
	}
#endif // MPI
	
	if (p1.lastcollision==t || p2.lastcollision==t) return;

	if (p1.number < N_tree_fixed) {
	  particles[c.p2].x += 2*boxsize;
	  fprintf(stderr,"Planet Collision\n");
	  return;
	}
	if (p2.number < N_tree_fixed) {
	  particles[c.p1].x += 2*boxsize;
	  fprintf(stderr,"Planet Collision\n");
	  return;
	}
	
	// Calculate new size distributions
	
	// Constants
	float fKE = 1.0;		// Fraction of kinetic energy used
	int rho = 3.0E3;		// Density (kg m^-3)
	float n0 = 2.8;			// Index of fragment distribution
	//float eqzero = 1.0e-20;	// Cutoff for "empty" bin
	double mconv = 1.98892e30; // conversion from solar units to kg
	float dconv = 1.5e11;	  // conversion from AU to m
	float vconv = 29866.;	  // conversion to m/s
	//float tconv = 0.159;      // conversion to years
	float vol1 = (4./3)*M_PI*pow(p1.r*dconv,3);	// Volume of swarm 1 (m^3)
	float vol2 = (4./3)*M_PI*pow(p2.r*dconv,3);  // Volume of swarm 2 (m^3)
	double Scon = 3.5e3;                    // strength regime coefficient
	double Gcon = 3.0e-8;                   // gravity regime coefficient
	float spow = -0.38;                    // strength regime exponent
	float gpow = 1.36;                     // gravity regime exponent
	float eta = -1.5;                      // supercatastrophic largest remnant exponent
	int extra = 30;
	int nowave = 1;
	
	// Build size bins
	double sizebins[numbins];
	for (int i=0; i<numbins; i++) {
	  sizebins[i] = pow(10,i*logbinsize-shift);
	}
	double bigbins[numbins+extra];
	for (int i=0; i<numbins+extra; i++) {
	  bigbins[i] = pow(10,i*logbinsize-shift-extra*logbinsize);
	}
	
	// Calculate minimum projectile size
	double mass[numbins];			// mass bins
	double bigvolume[numbins+extra];	// extrapolated volume bins
	double bigmass[numbins+extra];		// extrapolated mass bins
	double Qs[numbins+extra];		// internal binding energy
	double Qg[numbins+extra];		// gravitational binding energy
	double Qd[numbins+extra];		// minimum projectile kinetic energy
	double Ecol[numbins+extra][numbins+extra];  // collision energy
	double Qcol[numbins+extra][numbins+extra];  // collision energy/mass of target
	double Qsuper[numbins+extra][numbins+extra];// threshold for supercatastrophic fragmentation
	double Mlr[numbins+extra][numbins+extra]; // mass of largest fragment
	for (int i=0; i<numbins+extra; i++) {
		if (i>=extra) {
		  mass[i-extra] = rho*(4./3)*M_PI*pow(sizebins[i-extra]/2.,3); // small array of mass bins
		}
		bigvolume[i] = (4./3)*M_PI*pow(bigbins[i]/2.,3); // large array of volume bins
		bigmass[i] = rho*bigvolume[i]; // large array of mass bins
		Qs[i] = Scon*pow(bigbins[i]*100./2,spow); // calculate internal binding energy
		Qg[i] = Gcon*rho*pow(bigbins[i]*100./2,gpow); // calculate gravitational binding energy
		Qd[i] = (1./fKE)*(Qs[i]+Qg[i]); // calculate minimum projectile kinetic energy
	}
	
	// Velocity calculations
	struct ghostbox gb = c.gb;
	double x1 = p1.x + gb.shiftx; // x-coord of SP 1
	double y1 = p1.y + gb.shifty; // y-coord of SP 1
	double z1 = p1.z + gb.shiftz; // z-coord of SP 1
	double x2 = p2.x; // x-coord of SP 2
	double y2 = p2.y; // y-coord of SP 2
	double z2 = p2.z; // z-coord of SP 2
	double rp = p1.r+p2.r; // Sum of SP radii
	if (rp*rp < (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1)) return; // if the SPs are not overlapping, return
	double vx1 = p1.vx + gb.shiftvx; // x-vel of SP 1
	double vy1 = p1.vy + gb.shiftvy; // y-vel of SP 1
	double vz1 = p1.vz + gb.shiftvz; // z-vel of SP 1
	double vx2 = p2.vx; // x-vel of SP 2
	double vy2 = p2.vy; // y-vel of SP 2
	double vz2 = p2.vz; // z-vel of SP 2
	if ((vx2-vx1)*(x2-x1)+(vy2-vy1)*(y2-y1)+(vz2-vz1)*(z2-z1)>0) return; // if the SPs are not approaching each other, return
	
	// Calculate center-of-mass
	double m1 = p1.m; // mass of SP 1
	double m2 = p2.m; // mass of SP 2
	double cmx = (m1*vx1 + m2*vx2)/(m1 + m2); // x-coord of center of mass
	double cmy = (m1*vy1 + m2*vy2)/(m1 + m2); // y-coord of center of mass
	double cmz = (m1*vz1 + m2*vz2)/(m1 + m2); // z-coord of center of mass
	
	// Translate into center-of-mass frame
	double vxa = vx1 - cmx; // x-vel of SP 1 in CoM frame
	double vya = vy1 - cmy; // y-vel of SP 1 in CoM frame
	double vza = vz1 - cmz; // z-vel of SP 1 in CoM frame
	double vxb = vx2 - cmx; // x-vel of SP 2 in CoM frame
	double vyb = vy2 - cmy; // y-vel of SP 2 in CoM frame
	double vzb = vz2 - cmz; // z-vel of SP 2 in CoM frame
	double va = sqrt(vxa*vxa + vya*vya + vza*vza); // vel of SP 1 in CoM frame
	double vb = sqrt(vxb*vxb + vyb*vyb + vzb*vzb); // vel of SP 2 in CoM frame
	
	// Calculate pathlengths
	float relv = sqrt(pow(p1.vx-p2.vx,2)+pow(p1.vy-p2.vy,2)+pow(p1.vz-p2.vz,2));	// relative velocity of superparticles
	float pathlength1 = relv*(t-p1.lastcollision); // distance travelled by SP 1 since last collision 
	float pathlength2 = relv*(t-p2.lastcollision); // distance travelled by SP 2 since last collision 
	relv = relv*vconv; // relative velocity needs to be in m/s now
	float newrelv = relv; // relative velocity for each loop
	pathlength1 = pathlength1*dconv; // needs to be in m now
	pathlength2 = pathlength2*dconv; // needs to be in m now
	float pathlength;
	if (pathlength1 <= pathlength2) {
	  pathlength = pathlength1;
	} else {
	  pathlength = pathlength2;
	}

	// Correct pathlengths for superparticles from outside the disk
	if (pathlength1 > 5*mfp && collisions_Nlog > 10) {
	  pathlength1 = 5*mfp;
	  // fprintf(stderr,"Corrected pathlength1\n");
	}
	if (pathlength2 > 5*mfp && collisions_Nlog > 10) {
	  pathlength2 = 5*mfp;
	  //fprintf(stderr,"Corrected pathlength2\n");
	}

	// Initialize arrays for loss, survivors, and fragments
	double loss1[numbins];      // loss from SP 1
	double loss2[numbins];      // loss from SP 2
	double survivors1[numbins]; // survivors in SP 1
	double survivors2[numbins]; // survivors in SP 2
	double fragments1[numbins]; // fragments from SP 1
	double fragments2[numbins]; // frgments from SP 2
	// Initialize fragment arrays
	for (int i = 0; i < numbins; i++){
		fragments1[i] = 0.0;
		fragments2[i] = 0.0;
	}
	
	// Check for maximum optical depth
	double maxod = 0.0; // max optical depth
	double od1[numbins+extra]; // total optical depth for each bin in SP 1 due to SP 2
	double od2[numbins+extra]; // total optical depth for each bin in SP 2 due to SP 1
	double od1arr[numbins+extra][numbins+extra];	// [target][projectile]
	double od2arr[numbins+extra][numbins+extra];	// Prob of target in # hit by projectile in other
	double Eloss1 = 0.0; // total energy loss in SP 1
	double Eloss2 = 0.0; // total energy loss in SP 2
	double EL1[numbins],EL2[numbins]; // energy loss per bin in SP 1, SP 2
	double newdist1[numbins]; // new size distribution for SP 1
	double newdist2[numbins]; // new size distribution for SP 2
	// Start with new size distributions equal to initial size distributions
	for (int i=0; i<numbins; i++) {
		newdist1[i] = p1.sdist[i];
		newdist2[i] = p2.sdist[i];
	}

	float newpath1 = 0; // readjusted path length for SP 1
	float newpath2 = 0; // readjusted path length for SP 2
	float rempath1 = pathlength1; // remaining pathlength for SP 1
	float rempath2 = pathlength2; // remaining pathlength for SP 2
	int last = 0; // Flag for last loop
	
	//int numzeros1,numzeros2;
	//int numgood1,numgood2;
	//int k1,k2;
	struct line line1; // fit for extrapolated size distribution for SP 1
	struct line line2; // fit for extrapolated size distribution for SP 2
	double bigdist1[numbins+extra]; // extrapolated size distribution for SP 1
	double bigdist2[numbins+extra]; // extrapolated size distribution for SP 2

	double newm1,newm2; // new masses of SP 1 and SP 2
	double totm1,totm2; // new mass + lost dust mass for SP 1 and SP 2
	double KE; // total kinetic energy of SP 1 and SP 2
	double vc,vd,vxc,vyc,vzc,vxd,vyd,vzd; // new velocities for SP 1 and SP 2 in CoM frame
	float frac1,frac2; // fractional change in speed for SP 1 and SP 2

        int extrapolate1 = 0; // flag for whether to extrapolate SP 1
        int extrapolate2 = 0; // flag for whether to extrapolate SP 2
        int fitlength = numbins;    // number of bins to use for extrapolation
        double sdist1[fitlength];    // subset of SP 1 size dist to use for extrapolation
        double sdist2[fitlength];    // subset of SP 2 size dist to use for extrapolation
        double sizebins1[fitlength]; // subset of size bins to use for extrapolation
        double sizebins2[fitlength]; // subset of size bins to use for extrapolation

	int lr; // index of largest fragment
	double fragdistsum; // sum of fragment distribution (used to calculate B)
	double B; // normalization for fragment distribution
	double fragmentdist[numbins][numbins+extra][numbins];	// fragment distribution [target] [projectile] [fragments]
	
	double maxodsave;
	int empty1, empty2, empty;
	float eqzero = 1.0e-2;

	int numloop = 0; // number of loops
	do {  
	  numloop++; // track number of loops

	  ///////////////////////////////////////////////////////////////////////////////
	  // Extrapolate size distributions
	  // Check for empty bins
	  if (nowave) {
	    extrapolate1 = 1;
	    extrapolate2 = 1;
	    for (int i=0; i<fitlength; i++){
	      sdist1[i] = newdist1[i];
	      sdist2[i] = newdist2[i];
	      sizebins1[i] = sizebins[i];
	      sizebins2[i] = sizebins[i];
	      if (sdist1[i] <= 0) {
		extrapolate1 = 0; // If one of these bins is empty, don't extrapolate SP 1
	      }
	      if (sdist2[i] <= 0) {
		extrapolate2 = 0; // If one of these bins is empty, don't extrapolate SP 2
	      }
	    }
	  }
	  
	  // Extrapolate size distribution for SP 1 if necessary
	  if (extrapolate1) {
	    line1 = tools_linefit(sizebins1,sdist1,fitlength); // Fit to size distribution
	    for (int i=0; i<extra; i++) {
	      bigdist1[i] = pow(10.,line1.intercept)*pow(bigbins[i],line1.slope); // Build extrapolated distribution
	    }
	    // Error checking
	    if (isnan(bigdist1[0])) {
	      fprintf(stderr,"\nCollision between %.0f and %.0f yielded nan in %.0f -- bigdist\n",particles[c.p1].number,particles[c.p2].number,particles[c.p1].number);
	      for (int j=0; j<numbins; j++) {
		fprintf(stderr,"%e ",newdist1[j]);
	      }
	      fprintf(stderr,"\nslope = %f\n",line1.slope);
	      exit(0);
	    }
	  } else {
	    for (int i=0; i<extra; i++) {
	      bigdist1[i] = 0.0; // Don't build extrapolated distribution if not necessary
	    }
	  }
	  
	  
	  // Extrapolate size distribution for SP 2 if necessary
	  if (extrapolate2) {
	    line2 = tools_linefit(sizebins2,sdist2,fitlength); // Fit to size distribution
	    for (int i=0; i<extra; i++) {
	      bigdist2[i] = pow(10.,line2.intercept)*pow(bigbins[i],line2.slope); // Build extrapolated distribution
	    }
	    // Error checking
	    if (isnan(bigdist2[0])) {
	      fprintf(stderr,"\nCollision between %.0f and %.0f yielded nan in %.0f -- bigdist\n",particles[c.p1].number,particles[c.p2].number,particles[c.p2].number);
	      for (int j=0; j<numbins; j++) {
		fprintf(stderr,"%e ",newdist2[j]);
	      }
	      fprintf(stderr,"\nslope = %f\n",line2.slope);
	      exit(0);
	    }
	    
	  } else {
	    for (int i=0; i<extra; i++) {
	      bigdist2[i] = 0.0; // Don't build extrapolated distribution if not necessary
	    }
	  }
	  // Fill in the rest of the extrapolated distribution with the original distribution
	  for (int i=extra; i<numbins+extra; i++) {
	    bigdist1[i] = p1.sdist[i-extra]; 
	    bigdist2[i] = p2.sdist[i-extra];
	  }
	  	
	  ////////////////////////////////////////////////////////////////////////////
	  // Calculate optical depth
	  maxod = 0.0;
	  for (int i=0; i<numbins+extra; i++) {
	    od1[i] = 0.0;
	    od2[i] = 0.0;
	    for (int j=0; j<numbins+extra; j++) {
	      Ecol[i][j] = 0.5*bigmass[i]*bigmass[j]*newrelv*newrelv/(bigmass[i]+bigmass[j]);
	      Qcol[i][j] = Ecol[i][j]/bigmass[i];
	      Qsuper[i][j] = -2.0*Qd[i]*(0.1-1);
	      // Calculate size of largest remnant from collisional energy
	      if (Qcol[i][j] >= Qsuper[i][j]) {
		Mlr[i][j] = bigmass[i]*(0.1/pow(1.8,eta))*pow(Qcol[i][j]/Qd[i],eta);
	      } else {
		Mlr[i][j] = bigmass[i]*(-0.5*((Qcol[i][j]/Qd[i])-1)+0.5);
	      }
	      if (Mlr[i][j] < 0.0) {
		fprintf(stderr,"Error -- Mlr is less than zero \nEnding program\n");
		exit(0);
	      }
	      if (Mlr[i][j]/bigmass[i] < (1-pow(10,-3*logbinsize))) {
		od1arr[i][j] = bigdist2[j]*pow(bigbins[j]+bigbins[i],2)*(M_PI/4.)*rempath1/vol1;
		od2arr[i][j] = bigdist1[j]*pow(bigbins[j]+bigbins[i],2)*(M_PI/4.)*rempath2/vol2;
	      } else {
		od1arr[i][j] = 0.0;
		od2arr[i][j] = 0.0;
	      }
	      od1[i] += od1arr[i][j];
	      od2[i] += od2arr[i][j];
	    }
	    // Calculate max optical depth
	    if (od1[i] >= maxod && i >= extra) {
	      maxod = od1[i];
	    }
	    if (od2[i] >= maxod && i >= extra) {
	      maxod = od2[i];
	    }
	  }
	  if (numloop == 1) {
	    maxodsave = maxod;
	  }

	  // If max optical depth is greater than 1, only go that fraction of the pathlength
	  if (maxod > 1) {
	    newpath1 = 1.0*rempath1/maxod; // readjust pathlength
	    newpath2 = 1.0*rempath2/maxod; // readjust pathlength
	    // Readjust optical depths
	    for (int i=0; i<numbins+extra; i++) {
	      for (int j=0; j<numbins+extra; j++) {
		od1arr[i][j] = 1.0*od1arr[i][j]/maxod;
		od2arr[i][j] = 1.0*od2arr[i][j]/maxod;
	      }
	      od1[i] = 1.0*od1[i]/maxod;
	      od2[i] = 1.0*od2[i]/maxod;
	    }
	  // Otherwise, finish the loop
	  } else {
	    newpath1 = rempath1;
	    newpath2 = rempath2;
	    last = 1;
	  }

	  ///////////////////////////////////////////////////////////////////////////////////////
	  // Calculate fragment distributions
	  // Loop through each target
	  for (int i = 0; i < numbins; i++) {
	    // Loop through each projectile
	    for (int j = 0; j < numbins+extra; j++) {
	      //Mlr[i+extra][j] = 0.5*mass[i]; // Largest fragment is half the mass of the target
	      lr = 0;
	      // Find index of largest fragment
	      while (mass[lr] < Mlr[i+extra][j]) {
		lr ++;
	      }
	      if (lr > numbins) {
		fprintf(stderr,"Warning! Largest remnant outside array.\n");
		fprintf(stderr,"Mlr = %e\n",Mlr[i+extra][j]);
		fprintf(stderr,"bigmass = %e\n",bigmass[i+extra]);
		fprintf(stderr,"Diff = %e\n",Mlr[i+extra][j]-bigmass[i+extra]);
		exit(0);
	      }
	      // Normalize fragment distribution equal the mass of the target
	      fragdistsum = 0.0;
	      for (int k = 0; k <= (lr+extra); k++) {
		fragdistsum += pow(bigbins[k],3-n0);
	      }
	      B = 6.*mass[i]/(M_PI*rho*fragdistsum);
	      for (int k = 0; k < numbins; k++) {
		if (k <= lr) {
		  fragmentdist[i][j][k] = B*pow(sizebins[k],-1*n0);
		} else {
		  fragmentdist[i][j][k] = 0.0;
		}
	      }
	    }
	  }

	  for (int i=0; i<numbins; i++) {
	    fragments1[i] = 0.0;
	    fragments2[i] = 0.0;
	  }

	  /////////////////////////////////////////////////////////////////////////////
	  // Calculate survivors and fragments
	  for (int i=0; i<numbins; i++) {
	    EL1[i] = 0.0;
	    EL2[i] = 0.0;
	    // Calculate loss from each bin
	    loss1[i] = od1[i+extra]*newdist1[i];
	    loss2[i] = od2[i+extra]*newdist2[i];
	    for (int j=0; j<numbins+extra; j++) {
	      // Add to energy loss
	      // Each planetesimal in each collision loses half the collisional energy
	      EL1[i] += 0.5*Ecol[i][j]*loss1[i];
	      EL2[i] += 0.5*Ecol[i][j]*loss2[i];
	    }
	    // Total energy loss
	    Eloss1 += EL1[i];
	    Eloss2 += EL2[i];
	    // Number of survivors
	    survivors1[i] = newdist1[i] - loss1[i];
	    survivors2[i] = newdist2[i] - loss2[i];
	    // Add fragments to smaller bins
	    for (int k = 0; k < i; k++) {
	      for (int j = 0; j < numbins+extra; j++) {
		fragments1[k] += od1arr[i+extra][j]*newdist1[i]*fragmentdist[i][j][k];
		fragments2[k] += od2arr[i+extra][j]*newdist2[i]*fragmentdist[i][j][k];
	      }
	    }
	  }
	  // Calculate new size distributions 
	  for (int i=0; i<numbins; i++) {
	    newdist1[i] = survivors1[i]+fragments2[i];
	    newdist2[i] = survivors2[i]+fragments1[i];
	  }
	  // Calculate remaining pathlength
	  rempath1 = rempath1 - newpath1; 
	  rempath2 = rempath2 - newpath2;

	  // Calculate new total mass of each SP
	  newm1 = 0.0;
	  newm2 = 0.0;
	  for (int i=0; i<numbins; i++) {
	    newm1 += mass[i]*newdist1[i]/mconv; // now in solar masses
	    newm2 += mass[i]*newdist2[i]/mconv;
	  }

	  empty1 = 1;
	  for (int i=0; i<numbins; i++) {
	    if (newdist1[i] > eqzero) {
	      empty1 = 0;
	    } 
	  }
	  empty2 = 1;
	  for (int i=0; i<numbins; i++) {
	    if (newdist2[i] > eqzero) {
	      empty2 = 0;
	    } 
	  }
	  empty = 0;
	  if (empty1 || empty2) {
	    empty = 1;
	  }

	} while (! last && !empty && numloop < 101);  // While the optical depth is still greater than 1 and there is pathlength left, restart loop


	/////////////////////////////////////
	/////////////////////////////////////
	////////////END OF LOOP//////////////
	/////////////////////////////////////
	/////////////////////////////////////

	if (empty1) {
	  particles[c.p1].x += 2*boxsize;
	  fprintf(stderr,"Superparticle %d is empty\n",(int)particles[c.p1].number);
	}
	if (empty2) {
	  particles[c.p1].x += 2*boxsize;
	  fprintf(stderr,"Superparticle %d is empty\n",(int)particles[c.p2].number);
	}

	// Calculate new kinetic energy (in CoM frame for the two swarms)
	KE = 0.5*m1*va*va*mconv*vconv*vconv + 0.5*m2*vb*vb*mconv*vconv*vconv - (Eloss1+Eloss2); // Joules
	// Error checking
	if (KE < 0.0) {
	  fprintf(stderr,"\nERROR -- More energy used in collisions than available\n");
	  fprintf(stderr,"\nEloss1 = %e\n",Eloss1);
	  fprintf(stderr,"Eloss2 = %e\n",Eloss2);	
	  fprintf(stderr,"Einitial = %e\n",(0.5*m1*va*va+0.5*m2*vb*vb)*mconv*vconv*vconv);
	  fprintf(stderr,"KE = %e\n",KE);
	  exit(0);
	  KE = 0.0;
	}
	
	// Adjust velocities to conserve energy and momentum
	totm1 = (newm1/(newm1+newm2))*(m1+m2);
	totm2 = (newm2/(newm1+newm2))*(m1+m2);
	vc = sqrt(totm2*2*KE/(totm1*(totm1 + totm2)*mconv)); // m/s
	vd = (totm1/totm2)*vc; // m/s
	frac1 = va/vc;
	frac2 = vb/vd;
	if (frac1*frac2 == 0) {
	  vxc = vxa;
	  vyc = vya;
	  vzc = vza;
	  vxd = vxb;
	  vyd = vyb;
	  vzd = vzb;
	} else {
	  vxc = (vxa/frac1)/vconv; // now all in REBOUND velocity units
	  vyc = (vya/frac1)/vconv;
	  vzc = (vza/frac1)/vconv;
	  vxd = (vxb/frac2)/vconv;
	  vyd = (vyb/frac2)/vconv;
	  vzd = (vzb/frac2)/vconv;
	}
	
	// Calculate new relative velocity 
	newrelv = sqrt(pow(vxc-vxd,2)+pow(vyc-vyd,2)+pow(vzc-vzd,2));
	newrelv *= vconv;	  

	//struct orbit o1 = tools_p2orbit(particles[c.p1],particles[0]);
	//struct orbit o2 = tools_p2orbit(particles[c.p2],particles[0]);
	//fprintf(stderr,"%d loop(s): p1 = %e, p2 = %e, relv = %e, e1 = %e, e2 = %e, maxod = %e\n",numloop,pathlength1,pathlength2,relv,o1.e,o2.e,maxodsave);
	//if (numloop >= 10) {
	// fprintf(stderr,"maxod = %e\tpathlength1 = %e\tpathlength2=%e\n",maxodsave,pathlength1,pathlength2);
	//}
	//fprintf(stderr,"%d\t%e\t%e\t%e\n",numloop,pathlength1,pathlength2,maxodsave);

	if (isnan(newrelv) || numloop < 0) {
	  fprintf(stderr,"va = [%e,%e,%e]\n",vxa,vya,vza);
	  fprintf(stderr,"vb = [%e,%e,%e]\n",vxb,vyb,vzb);
	  fprintf(stderr,"pathlength1 = %e\n",pathlength1);
	  fprintf(stderr,"pathlength2 = %e\n",pathlength2);
	  fprintf(stderr,"mfp*7 = %e\n",mfp*7);
	  fprintf(stderr,"Collision between %d (%d) and %d (%d)\n",(int)particles[c.p1].number,(int)particles[c.p1].ncol,(int)particles[c.p2].number,(int)particles[c.p2].ncol);
	  fprintf(stderr,"Olddist1:\n");
	  for (int i=0; i<numbins; i++) {
	    fprintf(stderr,"%e,\t",p1.sdist[i]);
	  }
	  fprintf(stderr,"\n");
	  fprintf(stderr,"Olddist2:\n");
	  for (int i=0; i<numbins; i++) {
	    fprintf(stderr,"%e,\t",p2.sdist[i]);
	  }
	  fprintf(stderr,"\n");
	  fprintf(stderr,"Newdist1:\n");
	  for (int i=0; i<numbins; i++) {
	    fprintf(stderr,"%e\t",newdist1[i]);
	  }
	  fprintf(stderr,"\n");
	  fprintf(stderr,"Newdist2:\n");
	  for (int i=0; i<numbins; i++) {
	    fprintf(stderr,"%e\t",newdist2[i]);
	  }
	  fprintf(stderr,"\n");
	  fprintf(stderr,"newrelv = %e\n",newrelv);
	  fprintf(stderr,"numloop = %d\n",numloop);
	  fprintf(stderr,"maxod = %e\n",maxodsave);
	  exit(0);
	}

	// Translate back to original coordinate system
	double vx3 = vxc + cmx; // x-vel of SP 1 in REBOUND frame
	double vy3 = vyc + cmy; // y-vel of SP 1 in REBOUND frame
	double vz3 = vzc + cmz; // z-vel of SP 1 in REBOUND frame
	double vx4 = vxd + cmx; // x-vel of SP 2 in REBOUND frame
	double vy4 = vyd + cmy; // y-vel of SP 2 in REBOUND frame
	double vz4 = vzd + cmz; // z-vel of SP 2 in REBOUND frame

	// Record dust production
	//output_collisions_append("/discover/nobackup/enesvold/collisionsM3.0e0.4.txt",t,x1,y1,z1,vx3,vy3,vz3,totm1-newm1);
	//output_collisions_append("/discover/nobackup/enesvold/collisionsM3.0e0.4.txt",t,x2,y2,z2,vx4,vy4,vz4,totm2-newm2);	
	//output_collisions_append("/discover/nobackup/enesvold/relvelM1.0e0.3half.txt",t,x1,y1,z1,vx1,vy1,vz1,particles[c.p1].m);
	//output_collisions_append("/discover/nobackup/enesvold/relvelM1.0e0.3half.txt",t,x2,y2,z2,vx2,vy2,vz2,particles[c.p2].m);

	// Keep track of number of collisions and mean free path
	collisions_Nlog++;
	
	//fprintf(stderr,"%e: %f -- %f, %f -- %f\n",mfp,pathlength1/mfp,o1.e,pathlength2/mfp,o2.e);
	//fprintf(stderr,"%e, %e -- %e\n",pathlength1,pathlength2,mfp);
	// Don't count corrected pathlengths in mfp
	if (pathlength1 == 5*mfp) {
	  pathlength1 = mfp;
	}
	if (pathlength2 == 5*mfp) {
	  pathlength2 = mfp;
	}
	mfp = ((mfp*(collisions_Nlog-1)*2)+pathlength1+pathlength2)/(collisions_Nlog*2);
	//fprintf(stderr,"%e\n",mfp);
	//output_mfp("/discover/nobackup/enesvold/mfp.txt",t,mfp,pathlength1,pathlength2,numloop,p1.sdist[0],p2.sdist[0],maxodsave);
		
	// Apply the changes to the particles.
#ifdef MPI
	if (isloc==1){
#endif // MPI
	  particles[c.p2].vx = vx4;
	  particles[c.p2].vy = vy4;
	  particles[c.p2].vz = vz4;
	  particles[c.p2].m = newm2;
	  particles[c.p2].lastcollision = t;
	  for (int i=0; i<numbins; i++) {
	     	particles[c.p2].sdist[i] = newdist2[i];
	  }
	  particles[c.p2].ncol += 1;
#ifdef MPI
	}
#endif // MPI
	particles[c.p1].vx = vx3;
        particles[c.p1].vy = vy3; 
        particles[c.p1].vz = vz3;
	particles[c.p1].m = newm1;
	particles[c.p1].lastcollision = t; 
	for (int i=0; i<numbins; i++) {
		particles[c.p1].sdist[i] = newdist1[i];
	}
	particles[c.p1].ncol += 1;
	#endif // COLLISIONS_NONE
}

/**
 * @file 	tools.c
 * @brief 	Tools for creating distributions.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
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
#include "particle.h"
#include "main.h"
#include "tools.h"

double	tools_normaldistribution2_rsq;		/**< Used for speedup**/ 
double 	tools_normaldistribution2_v2;		/**< Used for speedup**/
int 	tools_normaldistribution2_ready = 0;	/**< Used for speedup**/

double tools_uniform(double min, double max){
	return ((double)rand())/((double)(RAND_MAX))*(max-min)+min;
}


double tools_powerlaw(double min, double max, double slope){
	double y = tools_uniform(0., 1.);
	return pow( (pow(max,slope+1.)-pow(min,slope+1.))*y+pow(min,slope+1.), 1./(slope+1.));
}

double tools_normal(double variance){
	if (tools_normaldistribution2_ready==1){
		tools_normaldistribution2_ready=0;
		return tools_normaldistribution2_v2*sqrt(-2.*log(tools_normaldistribution2_rsq)/tools_normaldistribution2_rsq*variance);
	}
	double v1,v2,rsq=1.;
	while(rsq>=1. || rsq<1.0e-12){
		v1=2.*((double)rand())/((double)(RAND_MAX))-1.0;
		v2=2.*((double)rand())/((double)(RAND_MAX))-1.0;
		rsq=v1*v1+v2*v2;
	}
	tools_normaldistribution2_ready = 1;
	tools_normaldistribution2_rsq   = rsq;
	tools_normaldistribution2_v2	  = v2;
	return 	v1*sqrt(-2.*log(rsq)/rsq*variance);
}

void tools_init_plummer(int _N, double mlow, double rfrac, int quiet, double scale, double* shift) {
	struct particle* _particles = calloc(_N,sizeof(struct particle));
	double scalefactor = (scale < 0 ?  16.0 / (3.0 * M_PI)  : scale);
	double inv_scalefactor = 1.0 / scalefactor;
	double sqrt_scalefactor = sqrt( scalefactor );
	rfrac *= scalefactor;          /* from VIRIAL to STRUCTURAL units */
	double mfrac = rfrac*rfrac*rfrac / pow(1.0 + rfrac*rfrac, 1.5);

	// Setup particles
	for (int i = 0; i < _N; i++) {
		_particles[i].m = 1.0/ (double) _N;
		double radius = 0.0;
		if (quiet==0){
			radius = 1.0 / sqrt( pow (tools_uniform(mlow,mfrac), -2.0/3.0) - 1.0);
		} else if (quiet==1) {
			double m_min = (i * mfrac)/(double)_N;
			double m_max = ((i+1) * mfrac)/(double)_N;
			radius = 1.0 / sqrt( pow (tools_uniform(m_min,m_max), -2.0/3.0) - 1.0);
		} else if (quiet==2) {
			double m_med = ((i+0.5) * mfrac)/(double)_N;
			radius = 1.0 / sqrt( pow (m_med, -2.0/3.0) - 1.0);
		} 	
		double theta = acos(tools_uniform(-1.0, 1.0));
		double phi = tools_uniform(0.0, 2.*M_PI);
		_particles[i].x = radius * sin( theta ) * cos( phi );
		_particles[i].y = radius * sin( theta ) * sin( phi );
		_particles[i].z = radius * cos( theta );
		double x = 0.0;
		double y = 0.1;
		while (y > x*x*pow( 1.0 - x*x, 3.5)) {
			x = tools_uniform(0.0,1.0);
			y = tools_uniform(0.0,0.1);
		}
		double velocity = x * sqrt(2.0) * pow( 1.0 + radius*radius, -0.25);
		theta = acos(tools_uniform(-1.0, 1.0));
		phi = tools_uniform(0.0,2.*M_PI);
		_particles[i].vx = velocity * sin( theta ) * cos( phi );
		_particles[i].vy = velocity * sin( theta ) * sin( phi );
		_particles[i].vz = velocity * cos( theta );

	}

	// Scale model and calculate center of mass.
	double w_x  = 0, w_y  = 0, w_z  = 0;
	double w_vx = 0, w_vy = 0, w_vz = 0;
	for (int i = 0; i < _N; i++) {
		_particles[i].x  *= inv_scalefactor;
		_particles[i].y  *= inv_scalefactor;
		_particles[i].z  *= inv_scalefactor;
		_particles[i].vx *= sqrt_scalefactor;
		_particles[i].vy *= sqrt_scalefactor;
		_particles[i].vz *= sqrt_scalefactor;
		w_x  += _particles[i].x;
		w_y  += _particles[i].y;
		w_z  += _particles[i].z;
		w_vx += _particles[i].vx;
		w_vy += _particles[i].vy;
		w_vz += _particles[i].vz;
	}
	double w_ins = 1./(double)_N;
	w_x *= w_ins;
	w_y *= w_ins;
	w_z *= w_ins;
	w_x -= shift[0];
	w_y -= shift[1];
	w_z -= shift[2];
	w_vx *= w_ins;
	w_vy *= w_ins;
	w_vz *= w_ins;
	w_vx -= shift[3];
	w_vy -= shift[4];
	w_vz -= shift[5];
	for (int i = 0; i < _N; i++) {
		_particles[i].x -= w_x;
		_particles[i].y -= w_y;
		_particles[i].z -= w_z;
		_particles[i].vx -= w_vx;
		_particles[i].vy -= w_vy;
		_particles[i].vz -= w_vz;
		particles_add(_particles[i]);
	}
	free(_particles);
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

#define TINY 1.0e-12
struct orbit tools_p2orbit(struct particle p, double cmass){
	struct orbit o;
	double h0,h1,h2,e0,e1,e2,n0,n1,n,er,vr,mu,ea;
	mu = G*(p.m+cmass);
	h0 = (p.y*p.vz - p.z*p.vy); 			//angular momentum vector
	h1 = (p.z*p.vx - p.x*p.vz);
	h2 = (p.x*p.vy - p.y*p.vx);
	o.h = sqrt ( h0*h0 + h1*h1 + h2*h2 );		// abs value of angular moment 
	double v = sqrt ( p.vx*p.vx + p.vy*p.vy + p.vz*p.vz );
	o.r = sqrt ( p.x*p.x + p.y*p.y + p.z*p.z );
	vr = (p.x*p.vx + p.y*p.vy + p.z*p.vz)/o.r;
	e0 = 1./mu*( (v*v-mu/o.r)*p.x - o.r*vr*p.vx );
	e1 = 1./mu*( (v*v-mu/o.r)*p.y - o.r*vr*p.vy );
	e2 = 1./mu*( (v*v-mu/o.r)*p.z - o.r*vr*p.vz );
 	o.e = sqrt( e0*e0 + e1*e1 + e2*e2 );		// eccentricity
	o.a = -mu/( v*v - 2.*mu/o.r );			// semi major axis
	o.P = 2.*M_PI*sqrt( o.a*o.a*o.a/mu );		// period
	o.inc = acos( h2/o.h ) ;				// inclination (wrt xy-plane)   -  Note if pi/2 < i < pi then the orbit is retrograde
	n0 = -h1;					// vector of nodes lies in xy plane => no z component
	n1 =  h0;		
	n = sqrt( n0*n0 + n1*n1 );
	er = p.x*e0 + p.y*e1 + p.z*e2;
	if (n<=1.e-30||o.inc<=1.e-30){			// we are in the xy plane
		o.Omega=0.;
		if (e1>=0.) { o.omega=acos(e0/o.e); }else{ o.omega = 2.*M_PI-acos(e0/o.e); }
	}else{
		if (e2>=0.) { o.omega=acos(( n0*e0 + n1*e1 )/(n*o.e)); }else{ o.omega=2.*M_PI-acos(( n0*e0 + n1*e1 )/(n*o.e)); }// pericenter = 0 if pericenter = ascending node
		if (n1>=0.) { o.Omega = acos(n0/n); }else{  o.Omega=2.*M_PI-acos(n0/n);} 					// longitude of ascending node in xy plane, measured from x axis
	//	if (isnan(o.Omega)||isinf(o.Omega)) o.Omega=0.;
	}
	o.f = er/(o.e*o.r);
	ea = (1.-o.r/o.a)/o.e;
	if (o.f>1.||o.f<-1.){				// failsafe
		o.f = M_PI - M_PI * o.f;
		ea  = M_PI - M_PI * ea;
	}else{
		o.f = acos(o.f);			// true anomaly = 0 if planet at pericenter
		ea  = acos(ea);				// eccentric anomaly
	}
	
	if (vr<0.) { 
		o.f=2.*M_PI-o.f;	
		ea =2.*M_PI-ea;
	}
	o.l = ea -o.e*sin(ea)+o.omega;			// mean longitude
	if (o.e<=1.e-10){ 				//circular orbit
		o.omega=0.;
		o.f=0.; 				// f has no meaning
		o.l=0.;
	}

	return o;

	/*
	// Compute the angular momentum H, and thereby the inclination INC.
	double Omega, a, e, M, E=0, f=0, omega, inc; 	// orbital paramaers

	double u; int ialpha;			// internals
	
	double hx = p.y*p.vz - p.z*p.vy;
	double hy = p.z*p.vx - p.x*p.vz;
	double hz = p.x*p.vy - p.y*p.vx;
	double h2 = hx*hx + hy*hy +hz*hz;
	double h  = sqrt(h2);
	inc = acos(hz/h);

	// Compute longitude of ascending node CAPOM and the argument of latitude u.
	double fac = sqrt(hx*hx + hy*hy)/h;
	
	if(fac < TINY ){
	  	Omega = 0.;
	  	u = atan2(p.y,p.x);
	  	if(fabs(inc - M_PI) < 10.*TINY){
			u = -u;
	  	}
	}else{
	  	Omega = atan2(hx,-hy); 
	  	u = atan2 ( p.z/sin(inc) , p.x*cos(Omega) + p.y*sin(Omega));
	}

	while(Omega < 0.) Omega = Omega + 2.*M_PI;
	
	while(u < 0.) u = u + 2.*M_PI;
	

	//  Compute the radius R and velocity squared V2, and the dot product RDOTV, the energy per unit mass ENERGY .

	double r = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
	double v2 = p.vx*p.vx + p.vy*p.vy + p.vz*p.vz;
	double vdotr = p.x*p.vx + p.y*p.vy + p.z*p.vz;
	double energy = 0.5*v2 - G*cmass/r;

	//  Determine type of conic section and label it via IALPHA
	if(fabs(energy*r/(G*cmass)) < sqrt(TINY)){
		ialpha = 0;
	}else{
	   	if(energy < 0.) ialpha = -1;
	   	if(energy > 0.) ialpha = +1;
	}

	// Depending on the conic type, determine the remaining elements

	// ELLIPSE :
	if(ialpha == -1){
		a = -0.5*G*cmass/energy;
	  	fac = 1. - h2/(G*cmass*a);
		double cape, w;
		if (fac > TINY){
			e = sqrt ( fac );
             		double face =(a-r)/(a*e);

			//... Apr. 16/93 : watch for case where face is slightly outside unity
             		if ( face > 1.){
                		cape = 0.;
             		}else{
                		if ( face > -1.){
                   			cape = acos( face );
				}else{
                   			cape = M_PI;
				}
			}
			
            		if ( vdotr < 0. ) cape = 2.*M_PI - cape;
			double cw, sw;
	    		cw = (cos( cape) -e)/(1. - e*cos(cape));
	    		sw = sqrt(1. - e*e)*sin(cape)/(1. - e*cos(cape));
	    		w = atan2(sw,cw);
	    		while(w < 0.) w = w + 2.*M_PI;
	  	}else{
	    		e = 0.;
	    		w = u;
	    		cape = u;
		}
		f = w;
		E = cape;
	  	M = cape - e*sin (cape);
	  	omega = u - w;
	  	while(omega < 0.) omega = omega + 2.*M_PI;
	  	omega = omega - floor(omega/(2.*M_PI))*2.*M_PI;
	}
	// HYPERBOLA :
	if(ialpha == 1){
	  	a = 0.5*G*cmass/energy;
	  	fac = h2/(G*cmass*a);
		double w, capf;
          	if (fac > TINY){
 	    		e = sqrt ( 1. + fac );
	    		double tmpf = (a+r)/(a*e);
            		if (tmpf < 1.0){
              			 tmpf = 1.0;
			}
	    		capf = log(tmpf + sqrt(tmpf*tmpf -1.));
	    		if ( vdotr < 0. ) capf = - capf;
			double cw,sw;
	    		cw = (e - cosh(capf))/(e*cosh(capf) - 1. );
	    		sw = sqrt(e*e - 1.)*sinh(capf)/(e*cosh(capf) - 1. );
	    		w = atan2(sw,cw);
	    		if(w < 0.) w = w + 2.*M_PI;
	  	}else{
	// we only get here if a hyperbola is essentially a parabola so we calculate e and w accordingly to avoid singularities
	    		e = 1.;
	    		double tmpf = 0.5*h2/(G*cmass);
	    		w = acos(2.*tmpf/r -1.);
	    		if ( vdotr < 0.) w = 2.*M_PI - w;
	    		tmpf = (a+r)/(a*e);
	    		capf = log(tmpf + sqrt(tmpf*tmpf -1.));
	  	}

	  	M = e * sinh(capf) - capf;
	  	omega = u - w;
	  	if(omega < 0.) omega = omega + 2.*M_PI;
		omega = omega - floor(omega/(2.*M_PI))*2.*M_PI;
	}

	// PARABOLA : ( NOTE - in this case we use "a" to mean pericentric distance)

	if(ialpha == 0){
		double w;
		a =  0.5*h2/(G*cmass);
	  	e = 1.;
	  	w = acos(2.*a/r -1.);
	  	if ( vdotr < 0.) w = 2.*M_PI - w;
	  	double tmpf = tan(0.5 * w);
	  	M = tmpf* (1. + tmpf*tmpf/3.);
	  	omega = u - w;
	  	if(omega < 0.) omega = omega + 2.*M_PI;
	  	omega = omega - floor(omega/(2.*M_PI))*2.*M_PI; 	 
	}
	
	struct orbit o; 
	o.Omega 	= Omega;
	o.omega 	= omega;
	o.M		= M;
	o.r		= r;
	o.f		= f;
	o.E		= E;
	o.inc		= inc;
	o.h		= h;
	o.e		= e;
	o.a		= a;
	return o;
	*/
	
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

#define LAPLACE_EPS 1.0e-14

double laplace(double s, int i, int j, double a){
	double as, term, sum, factor1, factor2, factor3, factor4;
	int k,q, q0;
	
	as = a*a;
	
	if(j<=i)     /* compute first term in sum */
    {
		factor4 = 1.0;
		for(k=0; k<j; k++)
			factor4 *= (i - k);
		sum = factor4;
		q0=0;
    }
	else
    {
		q0 = (j + 1 - i) / 2;
		sum = 0.0;
		factor4 = 1.0;
    }
	
	/* compute factors for terms in sum */
	
	factor1 = s;
	factor2 = s + i;
	factor3 = i + 1.0;
	for(q=1;q<q0;q++)   /* no contribution for q = 0 */
    {
		factor1 *= s + q;
		factor2 *= s + i + q;
		factor3 *= i + 1.0 + q;
    }
	
	term = as * factor1 * factor2 / (factor3 * q);
	
	/* sum series */
	
	while(term*factor4 > LAPLACE_EPS)
    {
		factor4 = 1.0;
		for(k=0;k<j;k++)
			factor4 *= (2*q + i - k);
		sum += term * factor4;
		factor1 += 1.0;
		factor2 += 1.0;
		factor3 += 1.0;
		q++;
		term *= as * factor1 * factor2 / (factor3 * q);
    }
	
	/* fix coefficient */
	
	for(k=0;k<i;k++)
		sum *= (s + ((double) k))/(((double) k)+1.0);
	
	if(q0 <= 0)
		sum *= 2.0 * pow(a, ((double) i));
	else
		sum *= 2.0 * pow(a, ((double) 2*q0 + i - 2));
	
	return(sum);
}

double tools_planckF(double T, double nu) {
  double h = 6.62606957e-34;
  double c = 2.99792458e8;
  double kB = 1.3806488e-23;

  double B = (2*h*pow(nu,3)/pow(c,2))/(exp(h*nu/(kB*T))-1);
  return(B);  
}

double tools_planckWL(double T, double lambda) {
  double h = 6.62606957e-34;
  double c = 2.99792458e8;
  double kB = 1.3806488e-23;

  double B = (2*h*pow(c,2)/pow(lambda,5))/(exp(h*c/(lambda*kB*T))-1);
  return(B);
}

void tools_histogram(double* array, double* hist, int nbins, double min, double max, int numpoints) {
  for (int i=0; i<nbins; i++){
    hist[i] = 0.0;
  }
  float binsize = (max-min)/nbins;
  for (int i=0; i<numpoints; i++){
    hist[(int) floor((array[i]-min)/binsize)]++;
  }
}

void tools_normarr(double*array, int size, double value) {
  double max = 0;
  for (int i=0; i<size; i++){
    if (array[i] > max){
      max = array[i];
    } 
  }

  double n = value/max;
  for (int i=0; i<size; i++){
    array[i] *= n;
  }
}

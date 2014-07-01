/**
 * @file 	tools.h
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
#ifndef TOOLS_H
#define TOOLS_H
#include "particle.h"
/**
 * Struct representing a Keplerian orbit.
 */
struct orbit {
	double a;
	double r;	// Radial distance from central object
	double h;	// Angular momentum
	double P;	// Orbital period
	double l;
	double e;
	double inc;
	double Omega; 	// longitude of ascending node
	double omega; 	// argument of perihelion
	double f; 	// true anomaly
};

struct line {
	double slope;
	double intercept;
};

/**
 * Calculates a random variable in a given range.
 * @param min Minimum value.
 * @param max Maximum value.
 */
double tools_uniform(double min, double max);

/**
 * Calculates a random variable drawn form a powerlaw distribution.
 * @param min Minimum value.
 * @param max Maximum value.
 * @param slop Slope of powerlaw distribution.
 */
double tools_powerlaw(double min, double max, double slope);

/**
 * Calculate a random number with normal distribution.
 * Algorithm by D.E. Knut, 1997, The Art of Computer Programmin, Addison-Wesley. 
 * @param variance Variance of normal distribution.
 * @return Random number with normal distribution (mean 0). 
 */
double tools_normal(double variance);

/**
 * This function sets up a Plummer sphere.
 * @details This function is based on a routine from the NEMO package, P. Teuben (1995).
 * For details on the implementation see the Appendix of Aarseth, Henon and Wielen (1974). 
 * @param _N Number of particles in the plummer sphere.
 * @param mlow Lower mass fraction cutoff (can be 0).
 * @param rfrac Upper radius cutoff (the Plummer sphere is formally an inifitely large object).
 * @param quiet Noisyness of the model, 0=noise, 1=medium, 2=quiet.
 * @param scale Scales the final model before adding it to the simulation.
 * @param shift Shift the entire sphere in position and velocity space (6 values). 
 */

void tools_init_plummer(int _N, double mlow, double rfrac, int quiet, double scale, double* shift);

/**
 * This function calculates x,y,z,vx,vy,vz given a,e,i,omega,OMEGA,M
 */

struct particle tools_orbit2p(float a, float e, float i, float omega, float OMEGA, float M, float Ms, float Mp);

/**
 * This function calculated orbital elements for a given particle. The center of
 * mass is assumed to be at the origin.
 * @param p Particle.
 * @param cmass Mass of the central object.
 * @return Orbital parameters. 
 */
struct orbit tools_p2orbit(struct particle p, double cmass);

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
 * This function calculates the Laplace coefficients and Leverrier derivatives
 *      j
 *  j  d    i
 * a  ---  b (a)
 *      j   s
 *    da
 * @author: Marc Kuchner
 * @param s 
 * @param i
 * @param j
 * @param a
 * @return sum Laplace coefficient
 */
double laplace(double s, int i, int j, double a);
 
/**
 * This function calculates the Planck function for temperature T and 
 * frequency nu.
 * @author: Erika Nesvold
 * @param T temperature
 * @param nu frequency
 * @return B(T)
 */
double tools_planckF(double T, double nu);

/**
 * This function calculates the Planck function for temperature T and
 * wavelength lambda.
 * @author: Erika Nesvold
 * @param T temperature
 * @param lambda wavelength
 * @return B(T)
 */
double tools_planckWL(double T, double lambda);

/**
 * This function builds a histogram of an input array.
 * @author: Erika Nesvold
 * @param array pointer to input array
 * @param histogram pointer to histogram array
 * @param nbins number of bins for histogram
 * @param min minimum value in array
 * @param max maximum value in array
 * @param numpoints number of elements in array
 */
void tools_histogram(double* array, double* hist, int nbins, double min, double max, int numpoints);

/**
 * This function normalizes an array so the max is a given value.
 * @author: Erika Nesvold
 * @param array pointer to array
 * @param size number of elements in array
 * @param value desired max value of array
 */ 

void tools_normarr(double* array, int size, double value);
 
 #endif 	// TOOLS_H

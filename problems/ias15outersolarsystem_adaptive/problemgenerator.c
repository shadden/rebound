#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <gmp.h>
#include <time.h>
#include <getopt.h>
#include "../../src/particle.h"

int N = 0;
struct particle* particles;
const double G	= 0.01720209895*0.01720209895;
double tmax;

int input_get_int(int argc, char** argv, const char* argument, int _default);
void output_binary(char* filename);

int main(int argc, char* argv[]){
	particles = calloc(sizeof(struct particle),1000); // more then needed
	int testcase = input_get_int(argc,argv,"testcase",0);
	
	switch (testcase){
		case 0: // Kozai  89.8
		{
			double semia = 0.01;

			struct particle star; 
			star.m  = 1;
			star.x  = 0; star.y  = 0; star.z  = 0; 
			star.vx = 0; star.vy = 0; star.vz = 0;
			particles[N++] = star;

			
			// The planet 
			struct particle planet; 
			planet.m  = 0.01*star.m;
			planet.x  = semia; 
			planet.y  = 0; 
			planet.z  = 0; 
			planet.vx = 0; 
			planet.vy = sqrt(G*star.m/semia); 
			planet.vz = 0;
			particles[N++] = planet;
			
			// The perturber
			struct particle perturber; 
			perturber.m  = 0.1*star.m;
			perturber.x  = semia*10.; 
			perturber.y  = 0; 
			perturber.z  = 0; 
			double inc_perturber = 89.9; // 80.
			perturber.vx = 0; 
			perturber.vy = cos(inc_perturber/180.*M_PI)*sqrt(G*(star.m+perturber.m)/perturber.x); 
			perturber.vz = sin(inc_perturber/180.*M_PI)*sqrt(G*(star.m+perturber.m)/perturber.x); 
			particles[N++] = perturber;
			tmax	= 1e4*365.*sqrt(semia*semia*semia/star.m);
			break;
		}
		case 1: // Kozai  80
		{
			double semia = 0.01;

			struct particle star; 
			star.m  = 1;
			star.x  = 0; star.y  = 0; star.z  = 0; 
			star.vx = 0; star.vy = 0; star.vz = 0;
			particles[N++] = star;

			
			// The planet 
			struct particle planet; 
			planet.m  = 0.01*star.m;
			planet.x  = semia; 
			planet.y  = 0; 
			planet.z  = 0; 
			planet.vx = 0; 
			planet.vy = sqrt(G*star.m/semia); 
			planet.vz = 0;
			particles[N++] = planet;
			
			// The perturber
			struct particle perturber; 
			perturber.m  = 0.1*star.m;
			perturber.x  = semia*10.; 
			perturber.y  = 0; 
			perturber.z  = 0; 
			double inc_perturber =  80.;
			perturber.vx = 0; 
			perturber.vy = cos(inc_perturber/180.*M_PI)*sqrt(G*(star.m+perturber.m)/perturber.x); 
			perturber.vz = sin(inc_perturber/180.*M_PI)*sqrt(G*(star.m+perturber.m)/perturber.x); 
			particles[N++] = perturber;
			tmax	= 1e4*365.*sqrt(semia*semia*semia/star.m);
			break;
		}
		case 2:
			break;
	}


	output_binary("particles.bin");

	FILE* of;
	
	// create files for mercury
	of  = fopen("big.in","w"); 
	fprintf(of,")O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n");
	fprintf(of,") Lines beginning with `)' are ignored.\n");
	fprintf(of,")---------------------------------------------------------------------\n");
	fprintf(of," style (Cartesian, Asteroidal, Cometary) = Cartesian\n");
	fprintf(of," epoch (in days) = 0.0\n");
	fprintf(of,")---------------------------------------------------------------------\n");
	for (int i=1;i<N;i++){
		fprintf(of,"NAME%i m=%.16e \n",i,particles[i].m);
		fprintf(of," %.16e %.16e %.16e\n",particles[i].x,particles[i].y,particles[i].z);
		fprintf(of," %.16e %.16e %.16e\n",particles[i].vx,particles[i].vy,particles[i].vz);
		fprintf(of," 0. 0. 0.\n");
	}
	fclose(of);
	
	of = fopen("param.in","w"); 
	fprintf(of,")O+_06 Integration parameters  (WARNING: Do not delete this line!!)\n");
	fprintf(of,") Lines beginning with `)' are ignored.\n");
	fprintf(of,")---------------------------------------------------------------------\n");
	fprintf(of,") Important integration parameters:\n");
	fprintf(of,")---------------------------------------------------------------------\n");
	fprintf(of," algorithm (MVS, BS, BS2, RADAU, HYBRID etc) = SCHEME\n");
	fprintf(of," start time (days)=  0d0\n");
	fprintf(of," stop time (days) = %e\n",tmax);
	fprintf(of," output interval (days) = %e\n",tmax);
	fprintf(of," timestep (days) = STEPSIZE\n");
	fprintf(of," accuracy parameter = EPSILON\n");
	fprintf(of,")---------------------------------------------------------------------\n");
	fprintf(of,") Integration options:\n");
	fprintf(of,")---------------------------------------------------------------------\n");
	fprintf(of," stop integration after a close encounter = no\n");
	fprintf(of," allow collisions to occur = no\n");
	fprintf(of," include collisional fragmentation = no\n");
	fprintf(of," express time in days or years = days\n");
	fprintf(of," express time relative to integration start time = no\n");
	fprintf(of," output precision = high\n");
	fprintf(of," < not used at present >\n");
	fprintf(of," include relativity in integration= no\n");
	fprintf(of," include user-defined force = no\n");
	fprintf(of,")---------------------------------------------------------------------\n");
	fprintf(of,") These parameters do not need to be adjusted often:\n");
	fprintf(of,")---------------------------------------------------------------------\n");
	fprintf(of," ejection distance (AU)= 1000\n");
	fprintf(of," radius of central body (AU) = 0\n");
	fprintf(of," central mass (solar) = %e\n",particles[0].m);
	fprintf(of," central J2 = 0\n");
	fprintf(of," central J4 = 0\n");
	fprintf(of," central J6 = 0\n");
	fprintf(of," < not used at present >\n");
	fprintf(of," < not used at present >\n");
	fprintf(of," Hybrid integrator changeover (Hill radii) = 3.\n");
	fprintf(of," number of timesteps between data dumps = 50000000\n");
	fprintf(of," number of timesteps between periodic effects = 10000000\n");
	fclose(of);

}


char* input_get_argument(int argc, char** argv, const char* argument){
	opterr = 0;
	optind = 1;
  	while (1) {
      		struct option long_options[] = {
	  		{NULL, required_argument, 0, 'a'},
			{0,0,0,0}
		};

		long_options[0].name = argument;

      		/* getopt_long stores the option index here.   */
      		int option_index = 0;
		//				short options. format abc:d::
      		int c = getopt_long (argc, argv, "", long_options, &option_index);

      		/* Detect the end of the options.   */
      		if (c == -1) break;

      		switch (c)
		{
			case 'a':
				return optarg;
				break;
			default:
				break;
		}
  	}
	return NULL;
}
int input_get_int(int argc, char** argv, const char* argument, int _default){
	char* value = input_get_argument(argc,argv,argument);
	if (value){
		return atoi(value);
	}
	return _default;
}

void output_binary(char* filename){
	FILE* of = fopen(filename,"wb"); 
	fwrite(&N,sizeof(int),1,of);
	fwrite(&tmax,sizeof(double),1,of);
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		fwrite(&(p),sizeof(struct particle),1,of);
	}
	fclose(of);
}

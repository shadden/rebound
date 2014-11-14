/**
*  Read in an integration parameter file.
*  Parameter files should have the format
*	1:	NUMBER OF PLANETS
*	2:	INPUT TYPE: 1 - CARTESIAN	2- ORBEL2D
*	3:	TIME STEP
*	4:	INTEGRATION LENGTH
*	5:	NUMBER OF TIMESTEPS BETWEEN DUMPS
*	6:	CENTRAL STAR MASS
*	7:	BOX SIZE
**/
int parse_param_data(FILE *fi,int* npl,int *inputType, double* dt,double* tFinal, int* N_dt_timing, double* Mcen, double* boxSize);


/**
*  Read in planet cartesian coordiantes for Npl planets.
*  Coordinate files should have a line for each planet with the form:
*		mass x y z vx vy vz
**/
void parse_cartesian_data(FILE *fi,int Npl,double Mcen);


/**
*  Read in planet 2-D orbital element coordiantes for Npl planets.
*  Coordinate files should have a line for each planet with the form:
*		mass a e omega M
**/
void parse_orbel2D_data(FILE *fi,int Npl,double Mcen);
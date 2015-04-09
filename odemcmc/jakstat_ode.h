#ifndef _ODEFUNDEF
#define _ODEFUNDEF
#include <iostream>
#include <iomanip>
#include <nvector/nvector_serial.h>
#include <assert.h>
#include <sundials_direct.h>

#define N_SPECIES 14
#define N_PARAMS 4
#define MAX_CONV_FAIL 1000000
#define MAX_STEPS 1000000000
#define MAX_ERRFAILS 15
#define MIN_STEPSIZE 0.000000000000010000000000000000
#define MAX_STEPSIZE 100000000000000.00000

#define OBS_EPO 1
#define FX_EPO 1


//EPO transient
int nepot = 16;
double epot[] = {0,2,4,6,8,10,12,14,16,18,20,25,30,40,50,60};
double epox[] = {0.01713,0.14500,0.24420,0.76590,1.00000,0.86050,0.78290,0.57050,0.62170,0.33100,0.33880,0.31160,0.05062,0.02504,0.01163,0.01243};

//EPO double peak
int nepot2 = 18;
double epot2[] = {0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,40,50,60};
double epox2[] = {0.052,0.445,1.000,0.373,0.406,0.412,0.484,0.553,0.792,0.786,0.624,0.610,0.505,0.214,0.166,0.125,0.083,0.042,0.000};

//EPO for sustained
int nepot3 = 18;
double epot3[] = {0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,40,50,60};
double epox3[] = {0.00,0.15,0.25,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75};

double sigma1[] = {0.0740,0.0500,0.0660,0.0700,0.0650,0.0510,0.0530,0.0510,0.0400,0.0400,0.0480,0.0520,0.0540,0.0550,0.0440,0.0710};
double sigma2[] = {0.0840,0.0460,0.0380,0.0320,0.0330,0.0370,0.0390,0.0400,0.0300,0.0280,0.0300,0.0310,0.0320,0.0400,0.0460,0.0820};

double interpolate(double t,int epotype){
	switch(epotype){
	case 1:
		for(int i=1;i<nepot;i++){
			if(epot[i] > t){
				double x = epox[i-1] + (t-epot[i-1])* (epox[i]-epox[i-1])/(epot[i]-epot[i-1]);
				return x;
				}
			}
		break;
	case 2:
		for(int i=1;i<nepot2;i++){
			if(epot2[i] > t){
				double x = epox2[i-1] + (t-epot2[i-1])* (epox2[i]-epox2[i-1])/(epot2[i]-epot2[i-1]);
				return x;
				}
			}
		break;
	case 3:
		for(int i=1;i<nepot3;i++){
			if(epot3[i] > t){
				double x = epox3[i-1] + (t-epot3[i-1])* (epox3[i]-epox3[i-1])/(epot3[i]-epot3[i-1]);
				return x;
				}
			}
		break;
	}
	assert(false);
	}

double measfun(double *data_in, double *x_in, int ns, int nt){
	double s = 0.0;
	double *m1,*m2,d1,d2;
	
	m1 = new double[nt-1];
	m2 = new double[nt-1];
	
	double m1max = 0.0;
	double m2max = 0.0;


	for(int t=0;t<nt-1;t++){
		m1[t] = x_in[t*N_SPECIES+1]+2*x_in[t*N_SPECIES+2];
		m2[t] = m1[t] + x_in[t*N_SPECIES];
		m1max = ((m1max > m1[t]) ? m1max : m1[t]);
		m2max = ((m2max > m2[t]) ? m2max : m2[t]);
		}
	
	for(int t=0;t<nt-1;t++){
		d1 = data_in[t*ns];
		d2 = data_in[t*ns+1];
		s += (d1-m1[t]/m1max)*(d1-m1[t]/m1max)/(sigma1[t+1]*sigma1[t+1]);
		s += (d2-m2[t]/m2max)*(d2-m2[t]/m2max)/(sigma2[t+1]*sigma2[t+1]);
		}
	
	delete[] m1;
	delete[] m2;
	
	return s;
	}

double vnuc = 450.0;
double vcyt = 1400.0;
	
int odefun(double t, N_Vector x_in, N_Vector dx_in, void *f_data){
	double *p = (double*) f_data;
	double *x, *dx;
	double epo;

	x = NV_DATA_S(x_in);
	dx = NV_DATA_S(dx_in);
	epo = interpolate(t,OBS_EPO);
	
	dx[0] = -p[0]*x[0]*epo	+ 2*(vnuc/vcyt)*p[3]*x[13];
	dx[1] =  p[0]*x[0]*epo	- p[1]*x[1]*x[1];
	dx[2] = -p[2]*x[2]		+ 0.5*p[1]*x[1]*x[1];
	dx[3] =  (vcyt/vnuc)*p[2]*x[2] - p[3]*x[3];
	
	dx[4] = p[3]*(x[3]-x[4]);
	dx[5] = p[3]*(x[4]-x[5]);
	dx[6] = p[3]*(x[5]-x[6]);
	dx[7] = p[3]*(x[6]-x[7]);
	dx[8] = p[3]*(x[7]-x[8]);
	dx[9] = p[3]*(x[8]-x[9]);
	dx[10] = p[3]*(x[9]-x[10]);
	dx[11] = p[3]*(x[10]-x[11]);
	dx[12] = p[3]*(x[11]-x[12]);
	dx[13] = p[3]*(x[12]-x[13]);
	
	return 0;
};

int odefunfx(double t, N_Vector x_in, N_Vector dx_in, void *f_data){
	double *p = (double*) f_data;
	double *x, *dx;
	double epo;

	x = NV_DATA_S(x_in);
	dx = NV_DATA_S(dx_in);
	epo = interpolate(t,FX_EPO);
	
	dx[0] = -p[0]*x[0]*epo	+ 2*(vnuc/vcyt)*p[3]*x[13];
	dx[1] =  p[0]*x[0]*epo	- p[1]*x[1]*x[1];
	dx[2] = -p[2]*x[2]		+ 0.5*p[1]*x[1]*x[1];
	dx[3] =  (vcyt/vnuc)*p[2]*x[2] - p[3]*x[3];
	
	dx[4] = p[3]*(x[3]-x[4]);
	dx[5] = p[3]*(x[4]-x[5]);
	dx[6] = p[3]*(x[5]-x[6]);
	dx[7] = p[3]*(x[6]-x[7]);
	dx[8] = p[3]*(x[7]-x[8]);
	dx[9] = p[3]*(x[8]-x[9]);
	dx[10] = p[3]*(x[9]-x[10]);
	dx[11] = p[3]*(x[10]-x[11]);
	dx[12] = p[3]*(x[11]-x[12]);
	dx[13] = p[3]*(x[12]-x[13]);
	
	return 0;
};
#endif
//Probabilistic verification on ODE models using MCMC
//Copyright (C) 2014  Benjamin M. Gyori

//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.

//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.

//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "cvode_sim.h"

#define PSI_TYPE 1

typedef int (*odefunType)(double, N_Vector, N_Vector, void*);
	double fxfun(double* p, odefunType odef){
		double xinit[] = {2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
		double ts[] = {0,2,4,6,8,10,12,14,16,18,20,24,30,40,50,60};
		int nt = 16;
		int ns = 14;
		double* x_out = new double[ns*nt];
		cvode_sim(ts,nt,xinit,p,odef,x_out);
		
		double *s = new double[nt];
		for(int t=0;t<nt;t++){
			s[t] = 0;
			for(int i=3;i<ns;i++){
				s[t] = s[t] + x_out[t*ns+i];
				}
			}
		bool peak, lastval;
		double f;
		switch(PSI_TYPE){
		case 1:
			peak = false;
			lastval = false;
			if(s[11]>=1) peak = true;
			if(s[nt-1]<=0.5) lastval = true;
			f = (peak && lastval);
			break;
		case 2:
			peak = false;
			lastval = false;
			for(int t=0;t<nt;t++){
				if(s[t]>=1) peak = true;
				}
			if((s[nt-1]>=0.5)&&(s[nt-1]<=1)) lastval = true;
			f = (peak && lastval);
			break;
		case 3:
			peak = false;
			lastval = false;
			for(int t=0;t<nt;t++){
				if(s[t]>=1.5) peak = true;
				if(peak){
					if(s[t]<1.5) lastval=false;
					}
				}
			f = (peak && lastval);
			break;
		}

		delete[] x_out;
		delete[] s;
		return f;
		}
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


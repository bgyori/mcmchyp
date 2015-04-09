#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <iostream>
#include <cstring>

struct experiment {
	experiment(double* xinit_in, int ns_x, double* dataTimes_in, double* dataSet_in, int ns, int nt) : ns(ns), nt(nt) {
		xinit = new double[ns_x];
		dataTimes = new double[nt];
		dataSet = new double[nt*ns];
		memcpy(xinit,xinit_in,ns_x*sizeof(double));
		memcpy(dataTimes,dataTimes_in,nt*sizeof(double));
		memcpy(dataSet,dataSet_in,ns*nt*sizeof(double));
		}
	
	~experiment(){
		delete[] xinit;
		delete[] dataTimes;
		delete[] dataSet;
		}
	
		int ns;
		int nt;
		double* xinit;
		double* dataTimes;
		double* dataSet;
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


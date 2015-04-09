#ifndef UNIFORM_PRIOR_H
#define UNIFORM_PRIOR_H

struct uniform_prior {
	uniform_prior(double* llim_in, double* ulim_in, int np_in)
		: llim(llim_in), ulim(ulim_in), np(np_in){}	
	
	double operator()(double* p){
		for(int i=0;i<np;i++){
			if((p[i] < llim[i]) || (ulim[i] < p[i])){
				return 0;
				}
			}
		return 1;
		}
	private:
		int np;
		double* llim;
		double* ulim;
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


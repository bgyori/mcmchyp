#ifndef ODE_MODEL_H
#define ODE_MODEL_H

#include <iostream>

#ifdef MODEL_JAKSTAT
	#include "jakstat_ode.h"
#endif

#include "fxfun.h"

typedef int (*odefunType)(double, N_Vector, N_Vector, void*);
typedef double (*measfunType)(double*,double*,int,int);
typedef double (*fxType)(double*,odefunType);

struct odeModel {
	template<class random_type,class prior_type>
	odeModel(random_type& ran,prior_type& p0) {
		odef = odefun;
		measf = measfun;
		ns = N_SPECIES;
		np = N_PARAMS;
		p = new double[np];
		p[0] = 2;p[1] = 10;p[2] = 0.15; p[3] = 2.5;
		}
	~odeModel(){
		delete[] p;
		}
		
	double* getstate(){return p;}
	int getdim(){return np;}
	
	void setstate(double* pn){
		for(int i=0;i<np;i++) {
			p[i] = pn[i];
			}
		}
	
	double fx(){
		return fxfun(p,odefunfx);
		}
		

	odefunType odef;
	measfunType measf;

	int ns, np;

	private:
		double *p;
		void printvec(){
			std::cout << "(";
			for(int i=0;i<np;i++){
				std::cout << p[i] << " ";
				}
			std::cout << ")" << std::endl;
			}
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


#ifndef GAUSS_UPDATE_H
#define GAUSS_UPDATE_H

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <utility>

using namespace std;

template<class model_type,
         class hamiltonian_type,
         class acceptance_rate_type>
class gauss_update {

public:

  gauss_update(model_type& m, 
                     hamiltonian_type& e,
                     acceptance_rate_type& r,
					double* s)
    : model(m), energy(e), rate(r), sigma(s) {
		int d = model.getdim();
		xn = new double[d];
		first = true;
		}
	
	~gauss_update(){
		delete[] xn;
		}

	template<class random_type>
	unsigned int operator()(random_type& ran) { 
		double r1,r2,ri;
		int d = model.getdim();
		double* x = model.getstate();

            
		for(int i=0;i<d;i++){
            ri = ran.randn();
			xn[i] = x[i]+sigma[i]*ri;
			}

		double ratio = energy(xn,&llh_old,&llh_new,first);	// Calculate probability ratio

		first = false;
		double rt = rate(ratio);		// Calculate probability of step
		double r = ran.randu();		
		if(r < rt){							// Accept step
			model.setstate(xn);
			llh_old = llh_new;
			return 1;
			}
		else {	// Reject step
			return 0;
			}
	}

private:
  model_type&               model;
  hamiltonian_type&         energy;
  acceptance_rate_type&     rate;
  double *sigma;
  double *xn;
  double llh_old,llh_new;
  bool first;
  ofstream fh;
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


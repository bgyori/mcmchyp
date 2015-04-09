#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <iostream>
#include <cmath>

template<class model_type,class prior_type>
struct hamiltonian {
	hamiltonian(model_type& m, prior_type& prf, experiment& expD): model(m), expData(expD), priorfun(prf){
		x_out = new double[expData.nt*model.ns];
		}
	~hamiltonian(){
		delete[] x_out;
		}
  
	// Returns ratio of probabilities
	double operator()(double* new_state,double* llh_old_in, double* llh_new_out, bool first) {
		double prior_old, llh_old, prior_new, llh_new;
		int ndim = model.getdim();
		
		// Get old state
		double *state = model.getstate();
		// Evaluate old prior
		prior_old = priorfun(state);
		// Calculate old likelihood
		if(first){
			*llh_old_in = getllh(state);
		}
		llh_old = *llh_old_in;
		
		
		// Evaluate new prior
		prior_new = priorfun(new_state);
		if (prior_new == 0) {
			return 0;
			}
		// Calculate new likelihood
		llh_new = getllh(new_state);
		*llh_new_out = llh_new;
		
		// Return ratio
		double ratio = (prior_new/prior_old)*exp(-0.5*(llh_new - llh_old));
		//printf("Ratio: %lf\n",ratio);
		
		return ratio;
		}
  

  
private:
	experiment& expData;
	prior_type& priorfun;
	model_type& model;
	double *x_out;

	double getllh(double* p_in) const {
		// Get experiment initial condition
		double* x_in = expData.xinit;

		// Get experiment data times
		double* ts = expData.dataTimes;
		int nt = expData.nt;

		// Simulate model
		double* x_out = new double[model.ns*nt];
		cvode_sim(ts,nt,x_in,p_in,model.odef,x_out);
		
		// Calculate sum-squared difference
		double s = model.measf(expData.dataSet,x_out+model.ns,expData.ns,nt);
		//std::cout << "s: " << s << std::endl;
		
		delete[] x_out;
		return s;
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


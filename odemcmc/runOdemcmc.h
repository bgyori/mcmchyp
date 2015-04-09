#include "rng_mersenne.h"
#include "odeModel.h"
#include "gauss_update.h"
#include "hamiltonian.h"
#include "mh_rate.h"
#include "monte_carlo.h"
#include "measurement.h"
#include "uniform_prior.h"

typedef hamiltonian<odeModel,uniform_prior> hamiltonian_type;
typedef gauss_update<odeModel,hamiltonian_type,mh_rate> updater_type;
typedef monte_carlo<updater_type,odeModel,measurement> mc_type;


void runOdemcmc(int nRuns, int nSteps,int nRelaxSteps,int nStepInterval, int seed,  
		experiment& expd, double* llim, double* ulim, int np, double* sigma,
		double* fxMean,double* fxVar,std::vector<double>& fx,std::string fname){
				
	mh_rate						rate;
    rng_mersenne ran(seed);

	if(!fxMean){
		uniform_prior prf(llim,ulim,np);
		odeModel		model(ran,prf);
		hamiltonian_type	h(model,prf,expd);
		updater_type	updater(model,h,rate,sigma);
		measurement		data(1);
		mc_type			mc(updater,model,data,fname);
		
		mc(ran,nRelaxSteps);
		mc(ran,nSteps,nStepInterval);
		data.getAll(fx);
		}
	else {
		for (int i=0;i<nRuns;i++){
			uniform_prior prf(llim,ulim,np);
			odeModel		model(ran,prf);
			hamiltonian_type		h(model,prf,expd);
			updater_type	updater(model,h,rate,sigma);
			measurement		data;
			mc_type			mc(updater,model,data,fname);

			mc(ran,nRelaxSteps);
			mc(ran,nSteps,nStepInterval);
			fxMean[i] = data.getMean();
			fxVar[i] = data.getVar();
		}
	}
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


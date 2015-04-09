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

#include "experiment.h"
#include "runOdemcmc.h"
#include <ctime>
#include <sstream>
#include <string>


int main(int argc, char* argv[]){
	int id;
	if(argc < 2){
		id = 1;
		}
	else {
		id = atoi(argv[1]);
		}


	int ns = 14;
	int np = 4;

	double expdata_xinit[] = {2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0};
	double expdata_datatimes[] = {0,2,4,6,8,10,12,14,16,18,20,25,30,40,50,60};
	double expdata_dataset[] = {0.3315,0.9275,0.8645,0.7923,0.9635,0.7778,0.9279,0.7053,0.8162,0.6522,0.7553,0.5894,0.7680,0.5894,0.8416,0.6377,0.7680,0.6425,0.8010,0.6908,0.7832,0.6908,0.8086,0.7585,0.4888,0.8068,0.2782,0.9275,0.2553,0.9710};

	double llim[] = {0,0,0,0};
	double ulim[] = {5,30,1,5};

	int nRuns 			= 1;
	int nSteps 			= 1e5;
	int nRelaxSteps		= 1e4;
	int nStepInterval	= 1;
	bool getAllfx 		= true;
	double sigma[] 		= {0.02,0.5,0.01,0.02};

	std::vector<double> fx;
	double fxMean[nRuns];
	double fxVar[nRuns];

	experiment expd = experiment(expdata_xinit,ns,expdata_datatimes,expdata_dataset,2,16);

	srand48(time(NULL));
	int seed = ceil(1+(drand48()*1000000.0)+id);
	std::cout << "Seed: " << seed << std::endl;

	std::clock_t tbegin = std::clock();

	std::stringstream fnamefxss;
	fnamefxss << "odemcmc-jakstat-" << nRuns <<"-"<< nSteps <<"-"<< nRelaxSteps <<"-"<< nStepInterval << "-"<< id << ".txt";
	std::cout << "Output file: " << fnamefxss.str() << std::endl;

	#ifdef LOG_FILE
		std::stringstream fnamelogss;
		fnamelogss << "odemcmc-jakstat-log-" << nRuns <<"-"<< nSteps <<"-"<< nRelaxSteps <<"-"<< nStepInterval << "-"<<sigma<<"-"<< id << ".txt";
		std::cout << "Log file: " << fnamelogss.str() << std::endl;
	#else
		std::stringstream fnamelogss;
		fnamelogss << "";
	#endif

	std:cout << nSteps << " steps, " << nRelaxSteps << " burn-in steps, " << nStepInterval << " interval" << std::endl;


    if (getAllfx){
		 runOdemcmc(1,nSteps,nRelaxSteps,nStepInterval,seed,expd,llim,ulim,np,sigma,0,0,fx,fnamelogss.str());
		 } 
	else {
		runOdemcmc(nRuns,nSteps,nRelaxSteps,nStepInterval,seed,expd,llim,ulim,np,sigma,fxMean,fxVar,fx,fnamelogss.str());
		}

 	std::clock_t tend = std::clock();
	double elapsed_secs = double(tend - tbegin) / CLOCKS_PER_SEC;

	
	double m = 0;
	ofstream fxfh;
	fxfh.open(fnamefxss.str().c_str());
	int nRealSteps = nSteps/nStepInterval;
	for(int i=0;i<nRealSteps;i++){
		fxfh << fx[i] << std::endl;
		m += fx[i]/nRealSteps;
		}
	fxfh.close();
	
	std::cout << "Mean: " << m << std::endl;

	std::cout << "Time taken: " << elapsed_secs << " s" << std::endl;

	return 0;
	}





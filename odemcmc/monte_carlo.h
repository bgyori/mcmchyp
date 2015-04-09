#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H
#include <iostream>

template<class updater_type,class model_type,class data_type>
struct monte_carlo {

	monte_carlo(updater_type& u,model_type& m,data_type& d,std::string &fn) 
				: updater(u),model(m),data(d) { 
		#ifdef LOG_FILE
			fname = fn;
			fh.open(fname.c_str());
			fh << setprecision(5);
		#endif
		dataAcc = new data_type();
		}

	~monte_carlo(){
		#ifdef LOG_FILE
			fh.close();
		#endif
		}

	template<class random_type>
	void operator()(random_type& ran, unsigned int mcs,
					unsigned int interval) {

		unsigned int acc;
		double fx;
		for(int i=0;i<mcs;i++){
 			acc = updater(ran);
			dataAcc.accumulate(acc);
			if((i%interval)==0){
				fx = model.fx();
				data.accumulate(fx);
				}
			else {
				fx = -1;
				}
			#ifdef LOG_FILE
				double* x = model.getstate();
				for(int j=0;j<model.np;j++) fh << x[j] << "\t";
				fh << std::endl;
			#endif
			}
		std::cout << "Acceptance ratio: " << dataAcc.getMean() << std::endl;
		}
  
	template<class random_type>
	void operator()(random_type& ran, unsigned int mcs) {
		for(unsigned int i=0;i<mcs;i++) {
			updater(ran);
		
			#ifdef LOG_FILE
				double* x = model.getstate();	
				for(int j=0;j<model.np;j++) fh << x[j] << "\t";
				fh << std::endl;
			#endif
			}
		}
 
	private:
		updater_type&    updater;
		model_type&      model;
		data_type       dataAcc;
		data_type&		data;
		std::ofstream fh;
		std::string fname;
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


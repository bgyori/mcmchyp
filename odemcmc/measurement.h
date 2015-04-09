#ifndef MEASUREMENT_H
#define MEASUREMENT_H

#include "statistics.h"
#include <vector>

struct measurement {
	measurement() {saveData=0;}
	measurement(bool sd) {saveData = sd;}
	
	void accumulate(double x) {
		if(saveData)
			dataAll.push_back(x);
		data.accumulate(x);
		}

	void getAll(std::vector<double>& vIn){
		vIn.swap(dataAll);
		}

	double getMean(){
		return data.getMean();
		}

	double getVar(){
		return data.getVar();
		}

	private:
		statistics data;
		std::vector<double> dataAll;
		int saveData;
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


#ifndef STATISTICS_H
#define STATISTICS_H

#include <vector>

struct _statData {
  double sum1;
  double sum2;
  double num;
  double mean;
  double stddev;
  double stde;
};

struct statistics {

  typedef _statData value_type;

  statistics(){
    data.sum1 = 0; data.sum2 = 0; data.num = 0; data.mean = 0;
    data.stddev = 0; data.stde = 0;
	}
  
  void accumulate(double a) {
      data.sum1  += a; 
      data.sum2  += a*a; 
      data.num   += 1;
      double m1 = data.sum1/data.num;
      double m2 = data.sum2/data.num;
      data.mean   = m1;
      data.stddev = sqrt(m2-m1*m1);
      data.stde = data.stddev / sqrt(data.num-1);
  }
  
  double getMean() const {
	return data.mean;
	}
  
  double getVar() const {
	double s1 = data.sum1;
	double s2 = data.sum2;
	double n = data.num;
	return s2-s1*(s1/n);
  }

  double getSum() const {
	return data.sum1;
	}

  private:
    value_type     data;
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


#ifndef RNG_MERSENNE_H
#define RNG_MERSENNE_H

#include <cmath>
#include "dSFMT.h"

struct rng_mersenne {

    rng_mersenne(){
        rng_mersenne(12345);
        }

    rng_mersenne(int seed){
        dsfmt_init_gen_rand(&dsfmt, seed);
        pi2 = atan(1.0)*8.0;
        }
    
    double randu(){
        return dsfmt_genrand_close_open(&dsfmt);
        }

    double randn(){
        static int second;
        static double r1,r2;
        if(second == 0){
            r1 = randu(); r2 = randu();
            second = 1;
            return cos(pi2*r1)*sqrt(-2.0*log(r2));
            }
        else {
            second = 0;
            return sin(pi2*r1)*sqrt(-2.0*log(r2));
            }
        }

    private:
        dsfmt_t dsfmt;
        double pi2;
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


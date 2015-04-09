#ifndef CVODE_SIM_H
#define CVODE_SIM_H

#include <iostream>

#include "cvode.h"
#include "nvector/nvector_serial.h"
#include "time.h"
#include "string.h"
#include "stdlib.h"
#include "cvode_dense.h"

typedef int (*odefunType)(double, N_Vector, N_Vector, void*);

int cvode_sim (double* ts,int nt,double* x0_in,double* param,odefunType f, double* x_out) {
	double tol[2]={1e-6,1e-10};	// Set tolerances here
	int k=0;
	double t0 = 0;
	double *x0_tmp;
	
	/* Structures required for CVode */
	int			flag;
	N_Vector	x_tmp = NULL;
	void		*cvode_mem = NULL;
	double	tret;
	long nsteps;
  

	/* Copy initial condition to local structure */
	x0_tmp = (double*)malloc(N_SPECIES*sizeof(double));
	memcpy(x0_tmp,x0_in,N_SPECIES*sizeof(double));
	x_tmp = N_VMake_Serial( N_SPECIES, x0_tmp );

	/* Start up CVode */
	cvode_mem = CVodeCreate( CV_BDF, CV_NEWTON );

	/* Initialise CVode */
	if ( CVodeInit( cvode_mem, f, t0, x_tmp ) != CV_SUCCESS ) {
		N_VDestroy_Serial( x_tmp );
		std::cerr << "ERROR: Failed to initialise CVode";
	}

	/* Specify tolerances */
	if ( CVodeSStolerances( cvode_mem, tol[0], tol[1] ) != CV_SUCCESS ) {
		N_VDestroy_Serial( x_tmp );
		CVodeFree( &cvode_mem );
		std::cerr << "ERROR: Failed to set tolerances";
	}

	CVodeSetMaxStep( cvode_mem, MAX_STEPSIZE );
	CVodeSetMinStep( cvode_mem, MIN_STEPSIZE );
	CVodeSetMaxConvFails( cvode_mem, MAX_CONV_FAIL );
	CVodeSetMaxNumSteps( cvode_mem, MAX_STEPS );
	CVodeSetMaxErrTestFails( cvode_mem, MAX_ERRFAILS );
  
	/* Attach dense linear solver module */
	if ( CVDense( cvode_mem, N_SPECIES ) != CV_SUCCESS ) {
		  N_VDestroy_Serial(x_tmp);
		  CVodeFree( &cvode_mem );
		  std::cerr << "ERROR: Failed to attach linear solver module";
	}
  
	/*  We need to pass our parameters to the ODE file */
	if ( CVodeSetUserData( cvode_mem, param ) != CV_SUCCESS ) {
		std::cerr << "ERROR: Failed passing parameters and initial conditions";
	}
  
	
	memcpy(x_out,x0_in,N_SPECIES*sizeof(double));
	
	for(int i=0;i<nt-1;i++){
		flag = CVode(cvode_mem, ts[i+1], x_tmp, &tret, CV_NORMAL);
		if (flag < 0) 
			std::cerr << "ERROR: CVode step error";
		k++;
		memcpy(x_out+(i+1)*N_SPECIES, &NV_DATA_S(x_tmp)[0], sizeof(double)*N_SPECIES);
	}

	N_VDestroy_Serial( x_tmp );

	/* Free CVode memory */
	CVodeFree( &cvode_mem );
	free(x0_tmp);

	return k;
}

#endif

//Probabilistic verification on ODE models using MCMC
//Copyright (C) 2014  Benjamin M. Gyori

// This source file is derived from Prediction Uncertainty Analysis (puaMAT)
// file mexG.c by Joep Vanlier

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


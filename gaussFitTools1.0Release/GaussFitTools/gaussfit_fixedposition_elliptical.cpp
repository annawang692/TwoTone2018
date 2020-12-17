/***%
%GaussFitTools 
%Copyright (C) 2010 The Chancellor, Masters and Scholars of the University of Oxford.
%Authors: Oliver J Britton, Seamus J Holden.
%
%This program is free software; you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation; either version 2 of the License, or
%(at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License (below) for more details.
%
*/
///Gaussian fitting LEVMAR CODE
///as described http://www.ics.forth.gr/~lourakis/levmar/levmar.pdf
///and distributed http://www.ics.forth.gr/~lourakis/levmar/

 // Uncompiled C++ code to make a .mex file for use in MATLAB
 // Fits to a Continuous Elliptical Gaussian model fitting parameters - A, sigma_x, sigma_y, b, theta with Xo and Yo fixed by input


     #include "gaussfit.h"

     using namespace std;




///SETUP

  char *variant = "fixedposition_elliptical";
//IMPORTANT - this defines which gauss fitting variant you are using
//for clarity - set it to the suffix of this gaussfitmex_suffix file
//although if it shares everything upto output with another variant
//(eg. free may want to just have a define block at the top of
//the file that defines the extra options.
//variant is used to choose what functions are run to set constraints etc.
 char *model = "elliptical"; //specifies if model is a circular or elliptical gaussian for purposes of output (eg brightness definition)

  const int m = 7;
  //total number of parameters - A, sigma_x, sigma_y, b, Xo, Yo, theta (should be 5 for circular and 7 for elliptical gaussians in 2D)
  const int background_int = 4;
  // background is the Nth variable in the Above List - need this to select correct value for background guess
  //if using an information criterion to filter


 //==========Define Function Pointers to Model and Jacobian functions==================
//These specify the functions that will be passed to the lsq solver to calculate values for the model function and its analytic jacobian
  void(*fit_model)(double *p, double *x, int m, int n, void *data ) = elliptical_gauss_f;
  void(*fit_jacobian)(double *p, double *jac, int m, int n, void *data) = elliptical_gauss_jac;
  //edit here to change model used, gauss_f/_jac for continuous gaussian, elliptical_gauss_f/_jac for an elliptical gaussian
//====================================================================================

 //========DATA RETRIEVAL FROM MATLAB============
  const int numEqualityConstraints = 2;
  //number of linear equality constraints eg. sigma = 1.5 for fixed width
  //generally this will be the number of variables you want fixed

  const int numBounds = 2 * 4; //2 * number of variables you want to define bounds for, to account for upper and lower bounding, in this case sigmaX, sigmaY, Xo, Yo
                               //A and b are hardcoded as 0 to 65535 (16^4)

  int numInputs = 1 + numEqualityConstraints + numBounds + 3;
  //number of inputs expected from Matlab  -the data array, equality constraints, bounds, whether to use the guess supplied,
  //the initial guess and options for the fit

  const int numDataEntries = numEqualityConstraints + numBounds; //data[] is an array that passes all constants needed to the constraint functions
                                //this lets the program know how many entries will be in data
                                //usually should be equal to numEqualityConstraints + numBounds
 //==============================================




//Main program starts here
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{

   gaussfit_main (nlhs, plhs, nrhs, prhs, variant, model, m, background_int,
                  fit_model, fit_jacobian,
                  numEqualityConstraints, numBounds, numInputs, numDataEntries);
//Run the main program

       return;
}

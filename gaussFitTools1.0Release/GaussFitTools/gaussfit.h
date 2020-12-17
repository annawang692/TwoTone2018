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
%GNU General Public License for more details.
%
*/
#ifndef GAUSSFIT_H_INCLUDED
#define GAUSSFIT_H_INCLUDED

     #include <stdlib.h>
     #include <stdio.h>
     #include <iostream>
     #include <string>
     #include <vector>
     #include <fstream>
     #include <cmath>
     #include "mex.h"

     #include "levmar.h"
     #include "gaussfit_main.cpp"
     #include "fitfunctions.cpp"  //contains the functions to be fitted and their analytic jacobians
     #include "initialfit.cpp"   //calculates initial guess from crude analysis of the data
     #include "constraints.cpp"                //sets the constraints on the solution

#endif // GAUSSFIT_H_INCLUDED


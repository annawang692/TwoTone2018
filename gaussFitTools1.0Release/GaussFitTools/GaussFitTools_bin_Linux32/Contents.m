%
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

%gaussfit_free - fits to a continuous circular Gaussian with all parameters
%free.
%Syntax is:
% [brightness, fit parameters, lsq] = gaussfit_free(Image, lower bound on
% sigma, upper bound on sigma, lower bound on Xo, upper bound on Xo, lower
% bound on Yo, upper bound on Yo, use supplied initial guess (true or
% false), initial guess as 5 element vector - (A,sigma, b, Xo, Yo), options
% vector - (display output on console true or false, epsilon1, 2,3, max
% number of iterations))
%
%gaussfit_fixedsigma - fits to a continuous circular Gaussian with
%sigma fixed. Syntax is:
% [brightness, fit parameters, lsq] = gaussfit_fixedsigma(Image, fixed value of sigma, lower bound on
% sigma, upper bound on sigma, lower bound on Xo, upper bound on Xo, lower
% bound on Yo, upper bound on Yo, use supplied initial guess (true or
% false), initial guess as 5 element vector - (A,sigma, b, Xo, Yo), options
% vector - (display output on console true or false, epsilon1, 2,3, max
% number of iterations))
% Lower and upper bounds on sigma should be set to 0 and 65535 or another
% large number as are already constraining sigma to a single value.
%
%gaussfit_fixedposition - fits to a continuous circular Gaussian with
%position (Xo, Yo) fixed. Syntax is:
% [brightness, fit parameters, lsq] = gaussfit_fixedposition(Image, fixed value of Xo, fixed value of Yo, lower bound on
% sigma, upper bound on sigma, lower bound on Xo, upper bound on Xo, lower
% bound on Yo, upper bound on Yo, use supplied initial guess (true or
% false), initial guess as 5 element vector - (A,sigma, b, Xo, Yo), options
% vector - (display output on console true or false, epsilon1, 2,3, max
% number of iterations))
% Lower and upper bounds on Xo, Yo should be set to 0 and 65535 or another
% large number as are already constraining them to a single value.

%gaussfit_fixedsigmaposition - fits to a continuous circular Gaussian with
%both sigma and position fixed. Syntax is:
% [brightness, fit parameters, lsq] = gaussfit_fixedsigmaposition(Image,fixed value of sigma,
% fixed value of Xo, fixed value of Yo, lower bound on sigma, upper bound on sigma, 
% lower bound on Xo, upper bound on Xo, lower  bound on Yo, upper bound on Yo, 
% use supplied initial guess (true or false),
% initial guess as 5 element vector - (A,sigma, b, Xo, Yo), options
% vector - (display output on console true or false, epsilon1, 2,3, max
% number of iterations))
% Lower and upper bounds on sigma, Xo, Yo should be set to 0 and 65535 or another
% large number as are already constraining them to a single value.

%(For all, Image must be a SQUARE matrix.)




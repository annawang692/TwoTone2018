function [phot_count a normChi2 ] = fixedGaussFit( Image, point_pos, crop_radius, maxwidth,minwidth, initguess,useCPPfit)
%	 [phot_count a normChi2 ] = fixedGaussFit( Image, point_pos, crop_radius, maxwidth,minwidth,initguess,useCPPfit)
% Inputs: 
%	frame
%	point_pos
%	crop_radius
%	maxwidth
% fit single 2d gaussian with fixed x, y position.
%
% ONLY ALLOW SQUARE SUBIMAGES OTHERWISE RETURN 0
% 
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 Alpha, released 101115
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%

if ~exist('useCPPfit','var')
  useCPPfit = true;
end

%set the fitting convergence tolerances
verbose = false;
epsilon1 =  10^-7; %gradient on lsq -
epsilon2 =  10^-9; %gradient on fitParams - the most important
epsilon3 = 0;  %absoluteValueLSQ - problem dependent - usually ~10^3 - NEVER USE IT!
maxIter = 30; %how fast it is in the absence of signal
if useCPPfit == true
  optionVector = [verbose, epsilon1,epsilon2,epsilon3,maxIter];
else
  curvefitoptions = optimset( 'lsqcurvefit');
  curvefitoptions = optimset( curvefitoptions,'Jacobian' ,'on','Display', 'off',  'TolX', epsilon2, 'TolFun', epsilon1,'MaxPCGIter',1,'MaxIter',maxIter);
end


[sizey sizex] = size(Image);
X0=point_pos(1);
Y0=point_pos(2);

%round X0, Y0 to use as matrix locations
X0_int = round(X0); 
Y0_int = round(Y0);
crop_radius = round(crop_radius); %radius should already be an integer anyway

% setup the limits of the cropped image
xstart =  X0_int-crop_radius;
xfinish = X0_int+crop_radius;
ystart =  Y0_int-crop_radius;
yfinish = Y0_int+crop_radius;
% check if any of the limits are out of bounds - if so, skip that point
if (xstart<1) || (xstart > sizex) ||  (xfinish<1) || (xfinish > sizex) ...
	|| (ystart<1) || (ystart > sizey) ||  (yfinish<1) || (yfinish > sizey) 
	
	a = [0 0 0]; 
	normChi2 = NaN;
else
	
	%crop to a small area around the point
	fitIm = Image( ystart:yfinish, xstart:xfinish);
	[sizeyFit sizexFit] = size(fitIm);
	% set up the point location in the cropped image coords
	X_POSim = X0-xstart+1;
	Y_POSim = Y0-ystart+1;


	%if an intial guess is not supplied, calculate it
	if all(initguess==0)
		background = min(fitIm(:));
		amplitude = max(fitIm(:));
		widthguess = widthEstimate(fitIm)/2;
		initguess = [amplitude widthguess background ];
	end
	
	%make sure the point position is within the limits of the cropped image
	if ( X_POSim < 1) || ( X_POSim > sizexFit) || ( Y_POSim < 1) || ( Y_POSim > sizeyFit)
		a = [0 0 0];
	else
		%do the fit 
		if useCPPfit == true %use the faster C++ fitting library
		  [a, normChi2] = fixedPosGaussFit_cpp(fitIm, X_POSim, Y_POSim , initguess, maxwidth , minwidth,optionVector);
		else %use the matlab version
		  [a, normChi2] = fixedgaussfit_matlab(fitIm, X_POSim, Y_POSim , initguess, maxwidth , minwidth,curvefitoptions);
		end
	end
end

%display( EXITFLAG);
%total photon count I = 2*pi*I0ab
I0 = a(1);
stdev = a(2);
BG0 = a(3);

phot_count = 2*pi*stdev^2*I0;
%keyboard;
%----------------------------------------------------------------------------------------------
function widthEst = widthEstimate(m)
%function to better estimate width for initial guess
[sizey sizex] = size(m);
vx = sum(m);
vy = sum(m');

vx = vx.*(vx>0);
vy = vy.*(vy>0);

x = [1:sizex];
y = [1:sizey];

cx = sum(vx.*x)/sum(vx);
cy = sum(vy.*y)/sum(vy);

sx = sqrt(sum(vx.*(abs(x-cx).^2))/sum(vx));
sy = sqrt(sum(vy.*(abs(y-cy).^2))/sum(vy));

widthEst = 0.5*(sx + sy) ;


%----------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------
function [a, normChi2] = fixedPosGaussFit_cpp(inputIm, X_POS, Y_POS,initguess , maxwidth, minwidth,optionVector)
%	[a, normChi2] = fixedPosGaussFit_cpp(inputIm, X_POS, Y_POS,initguess , maxwidth, minwidth,optionVector)
%
% ******wrapper for OB gaussFit C++ functions written 090915******
%
% fits point spread function, 
% F = (a(1)*exp(    -((X-X0)).^2+(Y-Y0).^2) /(2*a(2)^2)   ) + a(3))
%	initguess = [amplitude widthguess background ];
% Uses fixed position - X_POS, Y_POS

% [brightness, fit parameters, lsq] = gaussfit_fixedposition(Image, fixed value of Xo, fixed value of Yo, lower bound on
% sigma, upper bound on sigma, lower bound on Xo, upper bound on Xo, lower
% bound on Yo, upper bound on Yo, use supplied initial guess (true or
% false), initial guess as 5 element  - (A,sigma, b, Xo, Yo), options
% vector - (display output on console true or false, epsilon1, 2,3, max
% number of iterations))

%	initguess = [amplitude widthguess background ];
% initGuess5Vector - (A,sigma, b, Xo, Yo)

% ******wrapper for OB gaussFit C++ functions written 090915******
%AMPSCALEFACTOR =max(inputIm(:))/100;
%110419 Bug fix to prevent zero signal spikes
AMPSCALEFACTOR =max(inputIm(:))/5;
if AMPSCALEFACTOR <= 0 
  AMPSCALEFACTOR = 1;
end

%rescale the variables - to make magnitude of amplitude, background and width similar
inputIm = inputIm./AMPSCALEFACTOR;

useInitGuess = true;
initguess(1) = initguess(1)./AMPSCALEFACTOR;
initguess(3) = initguess(3)./AMPSCALEFACTOR;
initGuess5Vector = [initguess, X_POS, Y_POS];

[brightness, fitParams, lsq] = ...         
         gaussfit_fixedposition(inputIm,X_POS, Y_POS, minwidth, maxwidth,-65536, 65535, -65535, 65535, useInitGuess, initGuess5Vector, optionVector);
%set position bounds to 0 to 65535 so that position is only
      		 %constrained by equality constraints and not by
      		 %upper/lower bounds on X_POS, Y_POS

% modify fitParams to fit with "a" output syntax from twotoneMain 
%fitParams = (A,sigma, b, Xo, Yo)
% a = [A ,sigma,b]
a = fitParams(1:3);
if a(1)< 0 %We know that negative amplitude values are patently unphysical so ignore them
        a(1) = 0;
end

a(1) = a(1).*AMPSCALEFACTOR;
a(3) = a(3).*AMPSCALEFACTOR;
normChi2 = lsq; %this is the lsq normalised per pixel ie lsq_TOT.(Npix)^-1

%----------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------
function [a, normChi2] = fixedgaussfit_matlab(inputIm, X_POSim, Y_POSim , initguess, maxwidth , minwidth,curvefitoptions)

%set up the mesh, size of the input image for use in fitting
[sizey sizex] = size(inputIm);
num_pixels = sizey*sizex;
[X,Y]= meshgrid(1:sizex,1:sizey);
grid = [X Y];

%AMPSCALEFACTOR =max(inputIm(:))/100;
%110419 Bug fix to prevent zero signal spikes
AMPSCALEFACTOR =max(inputIm(:))/5;
if AMPSCALEFACTOR <= 0 
  AMPSCALEFACTOR = 1;
end
%rescale the variables - to make magnitude of amplitude, background and width similar
inputIm = inputIm./AMPSCALEFACTOR;
initguess(1) = initguess(1)./AMPSCALEFACTOR;
initguess(3) = initguess(3)./AMPSCALEFACTOR;


% Set fit limits on [amplitude width background]
lb = [0 minwidth 0];
ub = [65535 maxwidth 65535];
    
%do the fit
try
	[a, res] = ...
		lsqcurvefit(@(x, xdata) gauss2dw(x, xdata, X_POSim, Y_POSim), ...
		 initguess ,grid ,inputIm ,...
		  lb,ub,curvefitoptions);    
catch ME
	if strcmp(ME.identifier,'optim:snls:InvalidUserFunction') % supplied absolutely empty image!
		a = [0 0 0];
		res = 0;
	else
		rethrow(ME);
	end
end

if a(1)< 0 %We know that negative amplitude values are patently unphysical so ignore them
        a(1) = 0;
end

a(1) = a(1).*AMPSCALEFACTOR;
a(3) = a(3).*AMPSCALEFACTOR;
normChi2 = res/num_pixels; 	  


%-------------------------------------------------------------------------
%---------------------Fitting Subfunctions--------------------------------
%-------------------------------------------------------------------------

function [F J] = gauss2dw(a, data, X_POS, Y_POS)
% Used by the curve fitter to calculate values for a 2d gaussian
% with the x & y standard deviations equal
% and with fixed positions

%Initialise everything
[sizey sizex] = size(data);
sizex = sizex/2;

F = zeros(sizey, sizex);
X = F;
Y = F;

X = data(:,1:sizex);
Y = data(:,sizex+1:end);

% Only evaluate the exponential once:
expPart = exp( -((X-X_POS).^2+(Y-Y_POS).^2) /(2*a(2)^2) );

F = a(1)*expPart + a(3);

% compute the jacobian

% initialise everything
n = numel(F);
J = zeros(n,3); % initialise J
Ga1F = zeros(sizey, sizex);	
Ga2F = Ga1F;
Ga3F = Ga1F;

% Calculate the grad_a1(F),  grad_a2(F), grad_a3(F)
% see the printed notes on getGaussFit for derivation

Ga1F = expPart;

Ga2F = a(1)/a(2)^3 * ((X-X_POS).^2+(Y-Y_POS).^2) .* expPart;

Ga3F = ones(size(X));

% Form the jacobian, see the printed notes on getGaussFit for derivation
J = [Ga1F(:) Ga2F(:) Ga3F(:) ];



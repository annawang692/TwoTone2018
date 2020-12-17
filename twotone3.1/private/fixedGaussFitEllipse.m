function [phot_count a normChi2 ] = fixedGaussFitEllipse( im, point_pos, crop_radius, maxwidth,minwidth, initguess,useCPPfit)
%   [phot_count a normChi2 ] = fixedGaussFitEllipse( im, point_pos, crop_radius, maxwidth,minwidth,initguess,useCPPfit)
% Inputs: 
%  frame
%  point_pos
%  crop_radius
%  maxwidth
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
epsilon2 =  10^-9; %gradient on fitParam - the most important
epsilon3 = 0;  %absoluteValueLSQ - problem dependent - usually ~10^3 - NEVER USE IT!
maxIter  = 100; %how fast it is in the absence of signal
if useCPPfit == true
  optionVector = [verbose, epsilon1,epsilon2,epsilon3,maxIter];
else
  curvefitoptions = optimset( 'lsqcurvefit');
  curvefitoptions = optimset( curvefitoptions,'Jacobian' ,'on','Display', 'off',  'TolX', epsilon2, 'TolFun', epsilon1,'MaxPCGIter',1,'MaxIter',maxIter);
end





NUMFREEPARAMS = 7; % this is the number of free params in a, used for saving the fit params

[sizey sizex] = size(im);
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
  
  a =zeros(1,NUMFREEPARAMS); 
  normChi2 = NaN;
else
  
  %crop to a small area around the point
  fitIm = im( ystart:yfinish, xstart:xfinish);
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
      a =zeros(1,NUMFREEPARAMS); 
  else
    %do the fit 
    if useCPPfit == true %use the faster C++ fitting library
      %do the fit 
      [a, normChi2] =fixedPosEllipseGaussFit_cpp(fitIm, X_POSim, Y_POSim,initguess , maxwidth, minwidth,optionVector);
    else %use the matlab version
      %do the fit 
      [a, normChi2] =fixedPosEllipseGaussFit_matlab(fitIm, X_POSim, Y_POSim,initguess , maxwidth, minwidth,curvefitoptions);
    end
  end
end

%display( EXITFLAG);
%total photon count I = 2*pi*I0ab
% a: (A, sigma_x, sigma_y, b, Xo, Yo, theta )
I0 = a(1);
stdX = a(2);
stdY = a(3);
BG0 = a(4);

phot_count = 2*pi*stdX*stdY*I0;
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
function [fitParam, normChi2] = fixedPosEllipseGaussFit_cpp(im, X_POS, Y_POS,initguess , maxwidth, minwidth,optionVector)
%function [fitParam, normChi2] = fixedPosEllipseGaussFit_cpp(im, X_POS, Y_POS,initguess , maxwidth, minwidth,optionVector)
%
% ******wrapper for OB gaussFit C++ functions written 091215******
%
% fits point spread function, 
% F = (fitParam(1)*exp(    -((X-fitParam(4))).^2+(Y-fitParam(5)).^2) /(2*fitParam(2)^2)   ) + fitParam(3))
%  initguess = [amplitude widthguess background ];
%
% [brightness, fit parameters, lsq] = gaussfit_fixedposition(image, fixed value of Xo, fixed value of Yo, lower bound on
% sigma, upper bound on sigma, lower bound on Xo, upper bound on Xo, lower
% bound on Yo, upper bound on Yo, use supplied initial guess (true or
% false), initial guess as 5 element  - (A,sigma, b, Xo, Yo), options
% vector - (display output on console true or false, epsilon1, 2,3, max
% number of iterations))

%  initguess = [amplitude widthguess background ];
% initGuess7Vector - (A, sigma_x, sigma_y, b, Xo, Yo, theta )

sigmaXMin = minwidth;
sigmaXMax = maxwidth;
sigmaYMin = minwidth;
sigmaYMax = maxwidth;

%AMPSCALEFACTOR =max(im(:))/100;
%110419 Bug fix to prevent zero signal spikes
AMPSCALEFACTOR =max(im(:))/5;
if AMPSCALEFACTOR <= 0 
  AMPSCALEFACTOR = 1;
end
%AMPSCALEFACTOR =1;
%rescale the variables - to make magnitude of amplitude, background and width similar
im = im./AMPSCALEFACTOR;


initguess(1) = initguess(1)./AMPSCALEFACTOR; %amplitude
initguess(3) = initguess(3)./AMPSCALEFACTOR; %backgound
if ( (initguess(2) < minwidth) || (initguess(2) > maxwidth))%if given an out of bounds generate a
  initguess(2) = (maxwidth+minwidth)/2; %sensible one- careful this isnt too close to true val tho
end
thetaGuess = 0; 
initGuess7Vector = [initguess(1) initguess(2) initguess(2) initguess(3) X_POS Y_POS thetaGuess];
% A, sigma_x, sigma_y, b, Xo, Yo, theta 
if any(isnan(initGuess7Vector))
  useInitGuess = false;
else
  useInitGuess = true;
end


xMin = -65536;%even though its fixed need to define limits on position 
xMax = 65536;%just an idiosyncracy of the solver
yMin = -65536;
yMax = 65536;

[brightness,fitParam, lsq] =gaussfit_fixedposition_elliptical(im, X_POS, Y_POS,sigmaXMin,sigmaXMax,sigmaYMin,sigmaYMax,xMin,xMax,yMin,yMax , useInitGuess, initGuess7Vector, optionVector);

% modify fitParam to fit with "fitParam" output syntax from twotoneMain 
%fitParam = (A, sigma_x, sigma_y, b, Xo, Yo, theta )
% fitParam = [A ,sigma1, sigma2,b]
if fitParam(1)< 0 %We know that negative amplitude values are patently unphysical so ignore them
  fitParam(1) = 0;
end

fitParam(1) = fitParam(1).*AMPSCALEFACTOR;%amplitude
fitParam(4) = fitParam(4).*AMPSCALEFACTOR;%background
normChi2 = lsq; %this is the lsq normalised per pixel ie lsq_TOT.(Npix)^-1

%initGuess5Vector = [initguess(1) initguess(2)  initguess(3) X_POS Y_POS ]
%[brightness, fitParam, lsq] = ...         
%         gaussfit_fixedposition(im,X_POS, Y_POS, minwidth, maxwidth,xMin,xMax,yMin,yMax , useInitGuess, initGuess5Vector, optionVector);
%%set position bounds to 0 to 65535 so that position is only
%     %constrained by equality constraints and not by
%     %upper/lower bounds on X_POS, Y_POS
%
%% modify fitParam to fit with "a" output syntax from twotoneMain 
%%fitParam = (A,sigma, b, Xo, Yo)
%% a = [A ,sigma,b]
%a = fitParam(1:3);
%if a(1)< 0 %We know that negative amplitude values are patently unphysical so ignore them
%  a(1) = 0;
%end
%
%a(1) = a(1).*AMPSCALEFACTOR;
%a(3) = a(3).*AMPSCALEFACTOR;
%normChi2 = lsq; %this is the lsq normalised per pixel ie lsq_TOT.(Npix)^-1

%----------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------
function [a, normChi2] = fixedPosEllipseGaussFit_matlab(inputIm, X_POSim, Y_POSim , initguess, maxwidth , minwidth,curvefitoptions)
% fits point spread function, 
% F = (fitParam(1)*exp(    -(xprime).^2/(2*fitParam(2)^2)+(yprime).^2) /(2*fitParam(3)^2)   ) + fitParam(4))
%        
%           xprime = (X-Xo)*cos(theta) - (Y-Yo)*sin(theta);
%           yprime = (X-Xo)*sin(theta) + (Y-Yo)*cos(theta);
%
% extra fit params fitParam(7) = theta
% nb X0 and Y0 are saved as fitParam(5) and (6) even though they are fixed fact fixed
%
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
initguess(1) = initguess(1)./AMPSCALEFACTOR; %amplitude
initguess(3) = initguess(3)./AMPSCALEFACTOR; %backgound

%reformat initguess
thetaGuess = 0;
initGuess5Vector = [initguess(1) initguess(2) initguess(2) initguess(3) 0];

% Set fit limits on [amplitude widthx widthy background theta]
% dont set limits on theta but convert it to range 0->2pi afterwards
lb = [0 minwidth minwidth 0 -inf];
ub = [65535 maxwidth maxwidth 65535 inf];
    
%do the fit
try
	[fitParam, res] = ...
		lsqcurvefit(@(x, xdata) gauss2dw(x, xdata, X_POSim, Y_POSim), ...
		 initGuess5Vector ,grid ,inputIm ,...
		  lb,ub,curvefitoptions);    
catch ME
	if strcmp(ME.identifier,'optim:snls:InvalidUserFunction') % supplied absolutely empty image!
		fitParam = [0 0 0 0 0];
		res = 0;
	else
		rethrow(ME);
	end
end

a = [fitParam(1:4), X_POSim, Y_POSim, fitParam(5)];

if a(1)< 0 %We know that negative amplitude values are patently unphysical so ignore them
        a(1) = 0;
end

a(1) = a(1).*AMPSCALEFACTOR;
a(4) = a(4).*AMPSCALEFACTOR;
% convert theta to range 0->2 pi 
a(7) = mod(a(7),2*pi);
normChi2 = res/num_pixels; 	  

%-------------------------------------------------------------------------
%---------------------Fitting Subfunctions--------------------------------
%-------------------------------------------------------------------------

function [F J] = gauss2dw(a, data, X_POS, Y_POS)
% Used by the curve fitter to calculate values for a 2d gaussian
% with the x & y standard deviations equal
% and with fixed positions
% a(1) - A0
% a(2) - sX
% a(3) - sY
% a(4) - B
% a(5) - theta

%Initialise everything
[sizey sizex] = size(data);
sizex = sizex/2;

F = zeros(sizey, sizex);
X = F;
Y = F;

X = data(:,1:sizex);
Y = data(:,sizex+1:end);

xprime = (X-X_POS)*cos(a(5)) - (Y-Y_POS)*sin(a(5));
yprime = (X-X_POS)*sin(a(5)) + (Y-Y_POS)*cos(a(5));

% Only evaluate the exponential once:
expPart = exp( - ((xprime).^2 /(2*a(2)^2) + (yprime).^2 /(2*a(3)^2) ));

F = a(1)*expPart + a(4);

% compute the jacobian

% initialise everything
n = numel(F);
J = zeros(n,3); % initialise J
Ga1F = zeros(sizey, sizex);% dF/da(1)
Ga2F = Ga1F;% dF/da(2)
Ga3F = Ga1F;% dF/da(3)
Ga4F = Ga1F;% dF/da(4)
Ga5F = Ga1F;% dF/da(5)

% Calculate the grad_a1(F),  grad_a2(F), etc

Ga1F = expPart;

Ga2F = a(1).* expPart .*xprime.^2 .*a(2).^-3;% (A * e) * (pow(xprime,2) * pow(sigma_x,-3)); //dF/dsigma_x
Ga3F = a(1).* expPart .*yprime.^2 .*a(3).^-3;% (A * e) * (pow(yprime,2) * pow(sigma_y,-3)); //dF/dsigma_y

Ga4F = ones(size(X));
%dF/da(5) in c++ notation:
% (-A * e) *( (-xprime * pow(sigma_x,-2)) * ((X-Xo)*sin(theta) + (Y-Yo)*cos(theta)) + (yprime*pow(sigma_y,-2))*((X-Xo)*cos(theta) - (Y-Yo)*sin(theta)) );
Ga5F = -a(1).* expPart.* ...
      (	(-xprime).*a(2).^(-2).*((X-X_POS)*sin(a(5))+ (Y-Y_POS)*cos(a(5))) ... 
        + (yprime).*a(3).^(-2).*((X-X_POS)*cos(a(5))- (Y-Y_POS)*sin(a(5))) );


% Form the jacobian, see the printed notes on getGaussFit for derivation
J = [Ga1F(:) Ga2F(:) Ga3F(:) Ga4F(:) Ga5F(:)];




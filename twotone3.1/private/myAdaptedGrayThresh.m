function [level EM] = myAdaptedGrayThresh(image)
% function [level EM] = myAdaptedGrayThresh(image)
% returns gray thresh in correct range for integer type
% without failing for images packed into a small range
%
%assumes an image of double, uint16, uint8 etc.
% 
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 Alpha, released 101115
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%

maxIm = max(image(:));
minIm = min(image(:));
rangeIm = maxIm - minIm;

%class max & min:
imClass = class(image);
if strcmp(imClass,'double')
	maxClass = 1;
	minClass = 0;
else % its one of the integer types
	maxClass = intmax(imClass);
	minClass = intmin(imClass);
end
classRange = maxClass - minClass;


image = double(image);
dMinIm = double(minIm);
dMaxIm = double(maxIm);
%normalise the image
imageNormDouble = (image - dMinIm)/(dMaxIm - dMinIm);

[normThresh EM]= graythresh(imageNormDouble);

level = minIm + rangeIm*normThresh;

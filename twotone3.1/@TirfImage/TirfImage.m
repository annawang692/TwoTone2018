function TirfIm = TirfImage(impath, firstGreenFrame, imageLim, alternationPeriod)
% TirfImage class constructor.
% 	
%	TirfIm = TirfImage(impath, firstGreenFrame, imageLim, alternationPeriod)
%
% Description: 
%	Takes a tirf image splits it into red and green channels, and 
% performs an image transform on the green to map it to the red, using 
% transform, method for this is setTFORM_GrStack( TirfIm, TFORM)
%	TirfIm = TirfImage(impath) takes the tirf image and splits it
% but no transform
%
% Inputs:
% 	impath
%	firstGreenFrame
%	imageLim is a 2x4 matrix containing the limits of the green and red channels
%	  ie, [ greenXstart, greenXend, greenYstart, greenYend;
%		redXstart, redXend, redYstart, redYend];
%	alternationPeriod : if not supplied assume 2 if firstGreenFrame >0 ie 2colour alex
%			    if firstGreenFrame = 0 then cw so alternationPeriod ==1
%
% Data:	
%	GreenStack
%	RedStack
%	ImageLim   is a 2x4 matrix containing the limits of the green and red channels
%	  ie, [ greenXstart, greenXend, greenYstart, greenYend;
%		redXstart, redXend, redYstart, redYend];
%
%	FirstGreen - first green excitation frame (0/1/2) - 0 means
%			continuous excitiation
%	AlternationPeriod - period (in frames) of excitation.
%
% Public Methods:
%	- get[all]
%	- isGreenFrame(tirfim, frameno)      
%
% History:
% 	130208 	- First alpha complEte (SH)
%	180208 	- Function headers and help standardised (SH)
%		- Got rid of the redundant original image to save memory (SH)
%	200208  - Completely revised to make more memory efficient (SH)
%	240208	- Modified splitImage to automatically maximise the contrast (linearly)
%		  in the images. Very low contrast imput images were causing rounding errors
%		  which led to the stange results for background selection that
%		  ludo was getting (SH)
%	210708  - Takes the numFramesAvg as an opt. argument to only average the first n frames - THIS IS OBSOLETE 110808
% Author:
%	Seamus Holden
%
%
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 Alpha, released 101115
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%
%

if nargin == 0
	error('require > 0 arguments!')
end

if exist('alternationPeriod')
  TirfIm.AlternationPeriod = alternationPeriod;
elseif firstGreenFrame ==0 %ie cw excitation
  TirfIm.AlternationPeriod = 1;
else
 %default period is 2 ie 2colour alex
  TirfIm.AlternationPeriod = 2;
end

TirfIm.FirstGreen = firstGreenFrame;
TirfIm.ImageLim = imageLim;

%store the channels as ImageStack objects
TirfIm.GreenStack = ImageStack( impath, TirfIm.ImageLim(1,:));
TirfIm.RedStack   = ImageStack( impath, TirfIm.ImageLim(2,:));

%setup the class
TirfIm = class(TirfIm , 'TirfImage');


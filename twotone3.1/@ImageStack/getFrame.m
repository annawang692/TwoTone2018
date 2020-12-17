function Frame = getFrame(ImStack, n)
% 	function Frame = getFrame(ImStack, n)
% retrieve specified frame from ImStack
%
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 Alpha, released 101115
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only�? license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%

if n > ImStack.NumFrames || n < 1
	error('Attempted to access out of bounds frame');
end

if strcmp(ImStack.FileType,'fits')

%   if ImStack.NumFrames == 1 %ie if there is only one frame
% 	  startVector = double([ImStack.ImageLim(1), ImStack.ImageLim(3)]);
% 	  endVector =   double([ImStack.ImageLim(2), ImStack.ImageLim(4)]);
%   else
% 	  startVector = double([ImStack.ImageLim(1), ImStack.ImageLim(3), n]);
% 	  endVector =   double([ImStack.ImageLim(2), ImStack.ImageLim(4), n]);
%   end

  % the axis definitions in matlab for fits files are 
  % 90 degrees to the normal definition 
  % so rotate it
  % normal here means ImageJ
  %Frame = rot90(fits_read_image_subset(ImStack.ImPath, startVector, endVector));
  Frame = ImStack.movie(ImStack.ImageLim(3): ImStack.ImageLim(4),ImStack.ImageLim(1): ImStack.ImageLim(2),n);
elseif strcmp(ImStack.FileType,'tif')
  %y axis tif definition is inverted to normal definition (ijk versus xyz, so have to turn i-axis into y-axis)
  yMin = ImStack.NAXIS2 - ImStack.ImageLim(4) + 1;
  yMax = ImStack.NAXIS2 - ImStack.ImageLim(3) + 1;
  xMin = ImStack.ImageLim(1);
  xMax = ImStack.ImageLim(2);
  pixelRegion = {[yMin, yMax],[ xMin, xMax] };
  Frame = imread(ImStack.ImPath,n,'Info',ImStack.ImageInfo,'PixelRegion',pixelRegion);
else
  error('Unrecognised file type.');
end

assignin('caller', inputname(1), ImStack);



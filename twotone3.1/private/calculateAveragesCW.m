function [D,A, framesData] = calculateAveragesCW(tirfIm, avgFirst, avgLast)
%function [D,A, framesData] = calculateAveragesCW(tirfIm, avgFirst, avgLast)
% 
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 Alpha, released 101115
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%

green_input = getGreenStack(tirfIm);
red_input = getRedStack(tirfIm);

imtype = getImType(green_input);

numFrames = getNumFrames(green_input);
if (avgFirst > numFrames) 
  warning('twotone-calculateAveragesCW:outOfBoundsFrame', 'Error attempted to access out of bounds frame.\n Resetting avgFirst to 1\n');
  avgFirst =1;  
end
if (avgLast > numFrames)
  warning('twotone-calculateAveragesCW:outOfBoundsFrame', 'Error attempted to access out of bounds frame.\n Resetting avgLast to numFrames\n');
  avgLast = numFrames;
end

numGreenFrames	= avgLast - avgFirst + 1;
numRedFrames	= numGreenFrames;

if avgLast==0 || avgLast==0 || ...
    (avgFirst > avgLast) || (avgFirst > avgLast)
  error('Frame range is < 2 frames - not enough for an ALEX movie to be analysed');
end

% calculate the D average
D =  double(getFrame(green_input,avgFirst ))/numGreenFrames;
A = double( getFrame(red_input, avgFirst))/numRedFrames;
for i = avgFirst:avgLast
  greenimframe = double(getFrame(green_input, i));
  D = D + greenimframe/numGreenFrames;
  redimframe = double(getFrame(red_input, i));
  A = A + redimframe/numRedFrames;
end;
D = cast(D,imtype);
clear imframe;
clear green_input;
A = cast(A,imtype);
clear red_input;
clear imframe;

framesData.greenStart= avgFirst;
framesData.greenEnd = avgLast;
framesData.numGreenFrames = numGreenFrames;
framesData.redStart = avgFirst;
framesData.redEnd = avgLast;
framesData.numRedFrames = numRedFrames;


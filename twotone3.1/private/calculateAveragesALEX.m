function [DexDem,AexAem, DexAem, framesData] ...
  = calculateAveragesALEX(tirfIm, firstGreenFrame, avgFirst, avgLast)
% [DexDem,AexAem, DexAem] = calculateAveragesALEX(tirfIm, firstGreenFrame, avgFirst, avgLast)
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
  warning('twotone-calculateAveragesALEX:outOfBoundsFrame', 'Error attempted to access out of bounds frame.\n Resetting avgFirst to 1\n');
  avgFirst =1;  
end
if (avgLast > numFrames)
  warning('twotone-calculateAveragesALEX:outOfBoundsFrame', 'Error attempted to access out of bounds frame.\n Resetting avgLast to numFrames\n');
  avgLast = numFrames;
end

% work out which is the next green frame
%the next green frame, where first green frame i = 1,2
% and x is any integer is:
% nextGreen(x,i) = 2*round((x-i+1)/2) + i;
% the final green frame is
% finalGreen = nextGreen(x-1,i)
% similarly:
% nextRed(x,i) = nextGreen(x,3-i)
firstGreen = @(x,i) 2*round((x-i+2)/2) + i - 2;
finalGreen = @(x,i) firstGreen(x+1,i)-2;
firstRed = @(x,i) firstGreen(x,3-i);
finalRed = @(x,i) finalGreen(x,3-i);

%calculate the start and end frames for the specific values
greenStart = firstGreen(avgFirst,firstGreenFrame);
greenEnd = finalGreen(avgLast,firstGreenFrame);
redStart = firstRed(avgFirst,firstGreenFrame);
redEnd = finalRed(avgLast,firstGreenFrame);
if greenEnd==0 || redEnd==0 || ...
    (greenStart > greenEnd) || (redStart > redEnd)
  error('Frame range is < 2 frames - not enough for an ALEX movie to be analysed');
end

%calculate the number of frames which we are averaging over
numGreenFrames = numel(greenStart:2:greenEnd);
numRedFrames = numel(redStart:2:redEnd);

% calculate the DexDem average
DexDem =  double(getFrame(green_input,greenStart ))/numGreenFrames;
for i = greenStart:2:greenEnd
  imframe = double(getFrame(green_input, i));
  DexDem = DexDem + imframe/numGreenFrames;
end;
DexDem = cast(DexDem,imtype);
clear imframe;
clear green_input;

% calculate the AexAem and DexAem averages
% calculate DexAem
DexAem = double( getFrame(red_input, greenStart))/numGreenFrames;
for i = greenStart:2:greenEnd
  imframe = double(getFrame(red_input, i));
  DexAem = DexAem + imframe/numGreenFrames;
end;
DexAem = cast(DexAem,imtype);

% calculate AexAem
AexAem = double( getFrame(red_input, redStart))/numRedFrames;
for i = redStart:2:redEnd
  imframe = double(getFrame(red_input, i));
  AexAem = AexAem + imframe/numRedFrames;
end;
AexAem = cast(AexAem,imtype);
clear red_input;
clear imframe;

framesData.greenStart = greenStart;
framesData.greenEnd = greenEnd;
framesData.numGreenFrames = numGreenFrames;
framesData.redStart = redStart;
framesData.redEnd = redEnd;
framesData.numRedFrames = numRedFrames;


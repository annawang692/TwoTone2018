function returnVal = isGreenFrame(TirfIm,frameNo)
%function returnVal = isGreenFrame(TirfIm,frameNo)
%
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 Alpha, released 101115
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%

firstGreenFrame = TirfIm.FirstGreen;

if firstGreenFrame == 0
  firstGreenFrame =1;% treat CW as ALEX with FGF =1
end

if firstGreenFrame ==1
  if mod(frameNo,2)==1
    returnVal = true;
  else
    returnVal = false;
  end
else%for sencond frame green its the opposite way around
  if mod(frameNo,2)==0
    returnVal = true;
  else
    returnVal = false;
  end
end

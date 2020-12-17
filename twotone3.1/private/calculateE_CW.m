function fretData = calculateE_CW(frameNo,D,A,fretDataName)
%fretDataName should be: {'t','D','A','E'};
% 
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 Alpha, released 101115
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%

fretData = zeros(1, numel(fretDataName));

tRow = find(strcmp(fretDataName,'t'));
DRow = find(strcmp(fretDataName,'D'));
ARow = find(strcmp(fretDataName,'A'));
ERow = find(strcmp(fretDataName,'E'));

t = frameNo;
E = A/(D+A);
fretData(1,tRow) = t;
fretData(1,DRow)   = D;
fretData(1,ARow)   = A;
fretData(1,ERow)    = E;


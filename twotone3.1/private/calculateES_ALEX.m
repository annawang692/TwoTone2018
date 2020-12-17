function fretData = calculateES_ALEX(DexFrameNo,DD,DA,AD,AA, fretDataName)
%function fretData = calculateES_ALEX(DexFrameNo,DD,DA,AD,AA, fretDataName)
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

tdexRow = find(strcmp(fretDataName,'t_Dex'));
taexRow = find(strcmp(fretDataName,'t_Aex'));
DDRow = find(strcmp(fretDataName,'DD'));
DARow = find(strcmp(fretDataName,'DA'));
ADRow = find(strcmp(fretDataName,'AD'));
AARow = find(strcmp(fretDataName,'AA'));
ERow = find(strcmp(fretDataName,'E'));
SRow = find(strcmp(fretDataName,'S'));

t_Dex = DexFrameNo;
t_Aex = DexFrameNo + 1 ;
E = DA/(DD+DA);
S = (DD+DA)/(DD+DA+AA);
fretData(1,tdexRow) = t_Dex;
fretData(1,taexRow) = t_Aex;
fretData(1,DDRow)   = DD;
fretData(1,DARow)   = DA;
fretData(1,ADRow)   = AD;
fretData(1,AARow)   = AA;
fretData(1,ERow)    = E;
fretData(1,SRow)    = S;


function convertTwotoneToText(infilename,outfilename)
% convertTwotoneToText(infilename,outfilename)
%
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 , released 110426
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%
% FUNCTION: convertTwotoneToText
% DESCRIPTION:
%   Convert data to simple text file output.
% 
% INPUT: 
%   infilename  - Name of input data file
%   outfilename - Name to save output text file.
% OUTPUT: None

load(infilename,'twotoneData');
%header information
filename = twotoneData.results.movieInfo.name;
firstGreenFrame = twotoneData.results.movieInfo.firstGreenFrame;
nMolecule = numel(twotoneData.results.data);
if nMolecule >= 1
  nFrame = size(twotoneData.results.data(1).intensity,1);
else
  nFrame = NaN;
end
tformFile = twotoneData.settings.transformMatrixSettings.fileName{1};
fretDataName = twotoneData.settings.fretDataName;

fid = fopen(outfilename,'w');
fprintf(fid,'%% Twotone fitting results for FILE: %s\n',filename);
fprintf(fid,'%% Converted from mat-file MATFILE: %s\n',infilename);
fprintf(fid,'%% Using TFORMFILE: %s\n',tformFile);
fprintf(fid,'%% Number of frames NFRAME: %d\n',nFrame);
fprintf(fid,'%% Number of molecules: %d\n',nMolecule);
if firstGreenFrame == 0
  isALEX = 0;
  tRow = find(strcmp(fretDataName,'t'));
  DRow = find(strcmp(fretDataName,'D'));
  ARow = find(strcmp(fretDataName,'A'));
  ERow = find(strcmp(fretDataName,'E'));
  AdetDexRow = find(strcmp(twotoneData.settings.imageSettings.aDetChannelName,'D'));
  AdetAexRow = find(strcmp(twotoneData.settings.imageSettings.aDetChannelName,'A'));
else
  isALEX = 1;
  tdexRow = find(strcmp(fretDataName,'t_Dex'));
  taexRow = find(strcmp(fretDataName,'t_Aex'));
  DDRow = find(strcmp(fretDataName,'DD'));
  DARow = find(strcmp(fretDataName,'DA'));
  ADRow = find(strcmp(fretDataName,'AD'));
  AARow = find(strcmp(fretDataName,'AA'));
  ERow = find(strcmp(fretDataName,'E'));
  SRow = find(strcmp(fretDataName,'S'));
  AdetDDRow = find(strcmp(twotoneData.settings.imageSettings.aDetChannelName,'DexDem'));
  AdetAARow = find(strcmp(twotoneData.settings.imageSettings.aDetChannelName,'AexAem'));
  AdetDARow = find(strcmp(twotoneData.settings.imageSettings.aDetChannelName,'DexAem'));
end

aDetEccentricityRow = find(strcmp( twotoneData.results.aDetParamName,'eccentricity'));

fprintf(fid,'%% Movie is ALEX (1), CW (0), ISALEX: %d\n',isALEX);
fprintf(fid,'%% First green frame FIRSTGREEN: %d\n',firstGreenFrame);

%results
for i = 1:nMolecule
  fprintf(fid,'\n\n');
  fprintf(fid,'%% MOLECULE %d\n', i);
  nearestNeighbour = twotoneData.results.data(i).aDetData.nearestNeighborDist;
  fprintf(fid, '%% NEARESTNEIGHBOUR %f\n', nearestNeighbour);
  aDetPosData = twotoneData.results.data(i).aDetData.aDetPos;
  aDetEccentricity = twotoneData.results.data(i).aDetData.aDetParam(:,aDetEccentricityRow);
  fretData = twotoneData.results.data(i).fretData;

  if isALEX == 1
    t_Dex = fretData(:,tdexRow);
    t_Aex = fretData(:,taexRow);
    DD    = fretData(:,DDRow);  
    DA    = fretData(:,DARow);  
    AD    = fretData(:,ADRow);  
    AA    = fretData(:,AARow);  
    E     = fretData(:,ERow);   
    S     = fretData(:,SRow);   
    DD_XPOS =  aDetPosData(AdetDDRow,1);
    DD_YPOS =  aDetPosData(AdetDDRow,2);
    fprintf(fid,'%% Position (DexDem): XPOS %f, YPOS %f\n',DD_XPOS, DD_YPOS);
    fprintf(fid,'%% ADET_ECCENTRICITY_DD %f\n',aDetEccentricity(AdetDDRow));
    fprintf(fid,'%% ADET_ECCENTRICITY_DA %f\n',aDetEccentricity(AdetDARow));
    fprintf(fid,'%% ADET_ECCENTRICITY_AA %f\n',aDetEccentricity(AdetAARow));
    fprintf(fid,'%% T_DEX\tT_AEX\tDD\tDA\tAD\tAA\tE\tS\n');
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\n',...
      [t_Dex(:)';t_Aex(:)';DD(:)';DA(:)';AD(:)';AA(:)';E(:)';S(:)']);
  else
    t = fretData(:,tRow);
    D    = fretData(:,DRow);  
    A    = fretData(:,ARow);  
    E    = fretData(:,ERow);

    Dex_XPOS =  aDetPosData(AdetDexRow,1);
    Dex_YPOS =  aDetPosData(AdetDexRow,2);
    fprintf(fid,'%% Position (DexDem): XPOS %f, YPOS %f\n',Dex_XPOS, Dex_YPOS);
    fprintf(fid,'%% ADET_ECCENTRICITY_D %f\n',aDetEccentricity(AdetDexRow));
    fprintf(fid,'%% ADET_ECCENTRICITY_A %f\n',aDetEccentricity(AdetAexRow));

    fprintf(fid,'%% T\tD\tA\tE\n');
    fprintf(fid,'%d\t%d\t%d\t%f\n',...
      [t(:)';D(:)';A(:)';E(:)']);
  end
end
fprintf(fid,'\n');

fclose(fid);


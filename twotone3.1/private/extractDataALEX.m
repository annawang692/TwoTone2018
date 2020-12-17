function [twotoneDataOut] = extractDataALEX(fileName,twotoneData,clusteredData);
% function [logText] = extractDataALEX(fileName,twotoneData,clusteredData);
%   This function performs data analysis upon files
% 
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 Alpha, released 101115
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%


firstGreenFrame = twotoneData.results.movieInfo.firstGreenFrame;
imageLim = twotoneData.settings.imageSettings.channelImageLim;
alternationPeriod =twotoneData.settings.imageSettings.alternationPeriod;

tirfIm = TirfImage(fileName,firstGreenFrame, imageLim, alternationPeriod);
  
% run the analysis
if strcmp(twotoneData.settings.fitSettings.algorithm,'ring') %if its aperture photometry
  %TODO write aperturePhotometryALEX
  [twotoneData.results.data, twotoneData.results.fitParamName]= aperturePhotometryALEX( tirfIm,twotoneData, clusteredData);
else %profile fitting
  [twotoneData.results.data, twotoneData.results.fitParamName]= profileFittingPhotometryALEX( tirfIm,twotoneData, clusteredData);
end

twotoneDataOut = twotoneData;




function [data, fitParamName, fretDataName]= profileFittingPhotometryALEX( tirfim, twotoneData, clusteredData)
% function [data, fitParamName]= profileFittingPhotometryALEX( tirfim, twotoneData, clusteredData)
% INPUTS
%    % 2D circular Gaussian fit, fixed position
%    %twotoneData.settings.fitSettings.algorithm	= 'fixedGauss'; 
%      twotoneData.settings.fitSettings.param(1):	windowRadius 
%      twotoneData.settings.fitSettings.param(2):	minwidth
%      twotoneData.settings.fitSettings.param(3):	maxwidth
%    
%    % 2D elliptical Gaussian fit, fixed position
%    %twotoneData.settings.fitSettings.algorithm	= 'fixedGaussEllipse'; 
%     twotoneData.settings.fitSettings.param(1):	windowRadius 
%     twotoneData.settings.fitSettings.param(2):	minwidth
%     twotoneData.settings.fitSettings.param(3):	maxwidth
%    
%    % 2D elliptical Gaussian fit, free position
%    %twotoneData.settings.fitSettings.algorithm	= 'freeGaussEllipse'; 
%      twotoneData.settings.fitSettings.param(1):	windowRadius 
%      twotoneData.settings.fitSettings.param(2):	minwidth
%      twotoneData.settings.fitSettings.param(3):	maxwidth
%      twotoneData.settings.fitSettings.param(4):	poslim	%this is the number of pixels which the coordinate is allowed to move away from the autodetected position
%
% 
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 Alpha, released 101115
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only�? license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%

%for printing frame number to screen
LOOPVAL=5;

if exist('isGaussFitToolsInstalled','file') && isGaussFitToolsInstalled()==1
  useCPPfit = true;
  fitArg = {};
else
  useCPPfit = false;
  warning('GaussFitTools PSF fitting library not installed - expect extremely slow performance!');
  fitArg = {'useMatlabFit'};
end

if strcmp(twotoneData.settings.fitSettings.algorithm,'fixedGauss')
  windowRadius = twotoneData.settings.fitSettings.param(1);
  minwidth     = twotoneData.settings.fitSettings.param(2);
  maxwidth     = twotoneData.settings.fitSettings.param(3);
  NUMFREEPARAMS = 3; % this is the number of free params in a, used for saving the fit params
  fitParamName = {'A0', 's', 'BG'};%names of the free params
elseif strcmp(twotoneData.settings.fitSettings.algorithm,'fixedGaussEllipse')
  % fixed pos elliptical gaussian fitting
  windowRadius = twotoneData.settings.fitSettings.param(1);
  minwidth     = twotoneData.settings.fitSettings.param(2);
  maxwidth     = twotoneData.settings.fitSettings.param(3);
  NUMFREEPARAMS = 7; % this is the number of free params in a, used for saving the fit params
  fitParamName = {'A0', 'sx','sy', 'BG', 'X', 'Y', 'theta' };%names of the free params
elseif strcmp(twotoneData.settings.fitSettings.algorithm,'freeGaussEllipse')
  % free pos elliptical gaussian fitting
  windowRadius	= twotoneData.settings.fitSettings.param(1);
  minwidth     = twotoneData.settings.fitSettings.param(2);
  maxwidth	= twotoneData.settings.fitSettings.param(3);
  posLim	= twotoneData.settings.fitSettings.param(4);
  NUMFREEPARAMS = 7; % this is the number of free params in a, used for saving the fit params
  fitParamName = {'A0', 'sx','sy', 'BG', 'X', 'Y', 'theta' };%names of the free params
end

%a det param
aDetParamName = clusteredData(1).paramNames;
nAdetParam = numel(aDetParamName);

num_points = numel(clusteredData);
num_frames = getNumFrames(getRedStack(tirfim));

% is the image CW or ALEX
firstGreenFrame = getFirstGreen(tirfim);
% how many excitation channels do we have 
nFitChannel = twotoneData.settings.imageSettings.nChannel;
nADetChannel = twotoneData.settings.imageSettings.nADetChannel;
dFitCh = find(strcmp(twotoneData.settings.imageSettings.channelName ,'Dem'));
aFitCh = find(strcmp(twotoneData.settings.imageSettings.channelName ,'Aem'));
ddAdetCh = find(strcmp(twotoneData.settings.imageSettings.aDetChannelName ,'DexDem'));
aaAdetCh = find(strcmp(twotoneData.settings.imageSettings.aDetChannelName ,'AexAem'));
daAdetCh = find(strcmp(twotoneData.settings.imageSettings.aDetChannelName ,'DexAem'));

% initialise the intensities structure
aDetTempStruct = ...
  struct('aDetPos',zeros(nADetChannel,2),...
	  'aDetParam',zeros(nADetChannel,nAdetParam),...
	  'nearestNeighborDist', 0,...	    
	  'nParticle',zeros(nADetChannel,1));
dataTempStruct.aDetData = aDetTempStruct;
dataTempStruct.intensity = zeros(num_frames,nFitChannel);
dataTempStruct.fitParam = cellfun(@(x) zeros(num_frames,NUMFREEPARAMS), cell(nFitChannel,1),'UniformOutput',false);
nFretPairFrames = floor((double(num_frames) - double(~isGreenFrame(tirfim,1)) )/2);
fretDataName = twotoneData.settings.fretDataName;%should be  {t_Dex, t_Aex, DD, DA,AA, AA, E, S}
dataTempStruct.fretData = zeros(nFretPairFrames,numel(fretDataName));

data = repmat(dataTempStruct, num_points, 1);	

for i = 1:num_points
  %add all the autoDetection information
  data(i).aDetData.aDetPos(ddAdetCh,:)	 = clusteredData(i).dd;
  data(i).aDetData.aDetPos(aaAdetCh,:)	 = clusteredData(i).aa;
  data(i).aDetData.aDetPos(daAdetCh,:)	 = clusteredData(i).da;
  data(i).aDetData.nearestNeighborDist	 = clusteredData(i).NNcur;
  data(i).aDetData.nParticle(ddAdetCh)	 = clusteredData(i).ndd;
  data(i).aDetData.nParticle(aaAdetCh)	 = clusteredData(i).naa;
  data(i).aDetData.nParticle(daAdetCh)	 = clusteredData(i).nda;
  data(i).aDetData.aDetParam(ddAdetCh,:) = clusteredData(i).ddparams;
  data(i).aDetData.aDetParam(aaAdetCh,:) = clusteredData(i).aaparams;
  data(i).aDetData.aDetParam(daAdetCh,:) = clusteredData(i).daparams;
end

green_stack = getGreenStack(tirfim);
red_stack = getRedStack(tirfim);

% get the intensity plots for each point using this loop

% dont forget: intensities(i).point_pos is in (x,y) form whereas the stack matrices
% must be accessed in (i,j,k) form which is the opposite way round

%%
fprintf('Commencing data extraction\n');
kk =1;
for nn = 1:num_frames
  green_frame = double(getFrame( green_stack, nn));
  red_frame   = double(getFrame( red_stack, nn));

  for i = 1:num_points
    %green emission
    [data(i).intensity(nn,dFitCh) data(i).fitParam{dFitCh}(nn,:) ]=...
      gaussFit( green_frame,data(i).aDetData.aDetPos(ddAdetCh,:), twotoneData.settings.fitSettings,useCPPfit);
    %red emission
    if isGreenFrame(tirfim,nn)%use the DA position
    [data(i).intensity(nn,aFitCh) data(i).fitParam{aFitCh}(nn,:) ]=...
      gaussFit( red_frame,data(i).aDetData.aDetPos(daAdetCh,:), twotoneData.settings.fitSettings,useCPPfit);
    else %use the AA position
    [data(i).intensity(nn,aFitCh) data(i).fitParam{aFitCh}(nn,:) ]=...
      gaussFit( red_frame,data(i).aDetData.aDetPos(aaAdetCh,:), twotoneData.settings.fitSettings,useCPPfit);
    end

    %if this is a red frame > 1st frame, we can update the E,S calculation
    if ~isGreenFrame(tirfim,nn) && nn > 1
      DexFrameNo = nn - 1;
      DD = data(i).intensity(nn-1,dFitCh);
      DA = data(i).intensity(nn-1,aFitCh);
      AD = data(i).intensity(nn,dFitCh);
      AA = data(i).intensity(nn,aFitCh);
      data(i).fretData(kk,:) = calculateES_ALEX(DexFrameNo,DD,DA,AD,AA, fretDataName);
    end
  end

  if ~isGreenFrame(tirfim,nn) && nn > 1
    kk = kk + 1;
  end

  if nn==1 
    txtOut = sprintf('%d\t/\t%d',nn,num_frames);
    fprintf_r(txtOut,[],'reset');
  elseif nn==num_frames || mod(nn,LOOPVAL)==0 
    txtOut = sprintf('%d\t/\t%d',nn,num_frames);
    fprintf_r(txtOut,[]);
  end

end
fprintf('\nDone\n');

%----------------------------------------------------------------------------------------
function  [photCount  fitParamResult]  = gaussFit( im, pos, fitSettings,useCPPfit)
%function  [photCount  fitParamResult]  = gaussFit( im, pos, fitSettings,useCPPfit)
%run the fit, choosing between the appropriate methods. The fit functions are located in "gaussfitLib"

%DONT calculate an initial guess
initguess= 0;
if strcmp(fitSettings.algorithm,'fixedGauss')
  windowRadius = fitSettings.param(1);
  minwidth     = fitSettings.param(2);
  maxwidth     = fitSettings.param(3);
  [photCount  fitParamResult] = ...
    fixedGaussFit( im, pos, windowRadius, maxwidth,minwidth, initguess,useCPPfit);
elseif strcmp(fitSettings.algorithm,'fixedGaussEllipse')
  % fixed pos elliptical gaussian fitting
  windowRadius = fitSettings.param(1);
  minwidth     = fitSettings.param(2);
  maxwidth     = fitSettings.param(3);
  [photCount  fitParamResult] =...
    fixedGaussFitEllipse( im, pos, windowRadius, maxwidth,minwidth, initguess,useCPPfit);
elseif strcmp(fitSettings.algorithm,'freeGaussEllipse')
  % free pos elliptical gaussian fitting
  windowRadius = fitSettings.param(1);
  minwidth     = fitSettings.param(2);
  maxwidth     = fitSettings.param(3);
  posLim       = fitSettings.param(4);
  sigmaLim = [minwidth maxwidth];

  if useCPPfit == true 
    fitArg = {};
  else
    fitArg = {'useMatlabFit'};
  end

  [photCount  fitParamResult]=...
     freeGaussFitEllipse( im, pos, windowRadius,posLim, sigmaLim, initguess,fitArg{:});
else
  error('Fitting algorithm not recognised');
end


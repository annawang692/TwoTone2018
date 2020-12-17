function [data, fitParamName]= profileFittingPhotometryCW( tirfim, twotoneData, clusteredData)
% function [data, fitParamName]= profileFittingPhotometryCW( tirfim, twotoneData, clusteredData)
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
%      twotoneData.settings.fitSettings.param(2):	minwidth
%      twotoneData.settings.fitSettings.param(3):	maxwidth
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
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
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
  maxwidth     = twotoneData.settings.fitSettings.param(3);
  posLim	= twotoneData.settings.fitSettings.param(4);
  NUMFREEPARAMS = 7; % this is the number of free params in a, used for saving the fit params
  fitParamName = {'A0', 'sx','sy', 'BG', 'X', 'Y', 'theta' };%names of the free params
end

%a det param
aDetParamName = clusteredData(1).paramNames;
nAdetParam = numel(aDetParamName);

num_points = numel(clusteredData);
if exist('frames','var')
  frameRange = frames;
  num_frames = numel(frames);
else
  num_frames = getNumFrames(getRedStack(tirfim));
  frameRange = 1:num_frames;
end

% Image is CW
firstGreenFrame =0;

% how many excitation channels do we have 
nFitChannel = twotoneData.settings.imageSettings.nChannel;
nADetChannel = twotoneData.settings.imageSettings.nADetChannel;
dFitCh = find(strcmp(twotoneData.settings.imageSettings.channelName ,'D'));
aFitCh = find(strcmp(twotoneData.settings.imageSettings.channelName ,'A'));
dAdetCh = find(strcmp(twotoneData.settings.imageSettings.aDetChannelName ,'D'));
aAdetCh = find(strcmp(twotoneData.settings.imageSettings.aDetChannelName ,'A'));

% initialise the intensities structure
aDetTempStruct = ...
  struct('aDetPos',zeros(nADetChannel,2),...
	  'aDetParam',zeros(nADetChannel,nAdetParam),...
	  'nearestNeighborDist', 0,...	    
	  'nParticle',zeros(nADetChannel,1));
dataTempStruct.aDetData = aDetTempStruct;
dataTempStruct.intensity = zeros(num_frames,nFitChannel);
dataTempStruct.fitParam = cellfun(@(x) zeros(num_frames,NUMFREEPARAMS), cell(nFitChannel,1),'UniformOutput',false);
fretDataName = twotoneData.settings.fretDataName;%should be  {t_Dex, t_Aex, DD, DA,AA, AA, E, S}
dataTempStruct.fretData = zeros(num_frames,numel(fretDataName));

data = repmat(dataTempStruct, num_points, 1);	

for i = 1:num_points
  %add all the autoDetection information
  data(i).aDetData.aDetPos(dAdetCh,:) = clusteredData(i).d;
  data(i).aDetData.aDetPos(aAdetCh,:) = clusteredData(i).a;
  data(i).aDetData.nearestNeighborDist = clusteredData(i).NNcur;
  data(i).aDetData.nParticle(dAdetCh) = clusteredData(i).nd;
  data(i).aDetData.nParticle(aAdetCh) = clusteredData(i).na;
  data(i).aDetData.aDetParam(dAdetCh,:) = clusteredData(i).dparams;
  data(i).aDetData.aDetParam(aAdetCh,:) = clusteredData(i).aparams;
end

green_stack = getGreenStack(tirfim);
red_stack = getRedStack(tirfim);

% get the intensity plots for each point using this loop

% dont forget: intensities(i).point_pos is in (x,y) form whereas the stack matrices
% must be accessed in (i,j,k) form which is the opposite way round

%%
fprintf('Commencing data extraction\n');
lastFrame=max(frameRange); 
for nn = frameRange
  green_frame = double(getFrame( green_stack, nn));
  red_frame   = double(getFrame( red_stack, nn));

  for i = 1:num_points
    %green emission
    [data(i).intensity(nn,dFitCh) data(i).fitParam{dFitCh}(nn,:) ]=...
      gaussFit( green_frame,data(i).aDetData.aDetPos(dAdetCh,:), twotoneData.settings.fitSettings,useCPPfit);
    %red emission
    [data(i).intensity(nn,aFitCh) data(i).fitParam{aFitCh}(nn,:) ]=...
      gaussFit( red_frame,data(i).aDetData.aDetPos(aAdetCh,:), twotoneData.settings.fitSettings,useCPPfit);
    D = data(i).intensity(nn,dFitCh);
    A = data(i).intensity(nn,aFitCh);
    data(i).fretData(nn,:) = calculateE_CW(nn,D,A,fretDataName);
  end

  if nn==1 
    txtOut = sprintf('%d\t/\t%d',nn,lastFrame);
    fprintf_r(txtOut,[],'reset');
  elseif nn==lastFrame || mod(nn,LOOPVAL)==0 
    txtOut = sprintf('%d\t/\t%d',nn,lastFrame);
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


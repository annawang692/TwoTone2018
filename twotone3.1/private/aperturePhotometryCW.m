function [data, fitParamName]= aperturePhotometryCW( tirfim,  twotoneData,clusteredData)
%function [data, fitParamName]= aperturePhotometryCW( tirfim,  twotoneData,clusteredData)
%function intensities = aperturePhotometry( tirfim, clusteredData, innerCircleRadius, outerCircleRadius, ringTime)
% use the circle intensites methods to extract data
%
% convert every frame to a double to minimise roundoff error 140408 SH
%
%  ringFitParams = [circleAvg; ringAvg; totalCirclePixels; totalRingPixels];
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
LOOPVAL=20;
 
NUMFREEPARAMS = 4; % this is the number of free params in a, used for saving the fit params
fitParamName = {'circleAvg','ringAvg', 'totalCirclePixels', 'totalRingPixels'};%names of the free params
innerCircleRadius = twotoneData.settings.fitSettings.param(1);
outerCircleRadius = twotoneData.settings.fitSettings.param(2);
ringTime = twotoneData.settings.fitSettings.param(3);

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

% image is CW
firstGreenFrame = 0;

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
greenFrameSize = size(getFrame( green_stack, 1));
redFrameSize   = size(getFrame(   red_stack, 1));

for i = 1:num_points
  
  dPos = data(i).aDetData.aDetPos(dAdetCh,:);
  aPos = data(i).aDetData.aDetPos(aAdetCh,:);

  [DringIndex{i} DcircleIndex{i}] = ...
    calculateApertureIndex(greenFrameSize,dPos, innerCircleRadius, outerCircleRadius);

  [AringIndex{i} AcircleIndex{i}] = ...
    calculateApertureIndex(redFrameSize,aPos , innerCircleRadius, outerCircleRadius);
end

% get the intensity plots for each point using this loop

% dont forget: intensities(i).point_pos is in (x,y) form whereas the stack matrices
% must be accessed in (i,j,k) form which is the opposite way round
%
% if we are doing time averaged ring subtraction, initialize the totalBG matrices
% for speed
greenTotalBG = zeros(num_points, num_frames);
redTotalBG = zeros(num_points, num_frames);

fprintf('Commencing data extraction\n');
for nn = 1:num_frames

%  %timing test
%  tic;  
%  TOC = toc;
%  tic

  green_frame = double(getFrame( green_stack, nn));
  red_frame   = double(getFrame(   red_stack, nn));
  
  for i = 1:num_points
    
    subtract_bg = false; % no local ring subtraction
      
    %green emission
    [data(i).intensity(nn,dFitCh) greenTotalBG(i, nn) data(i).fitParam{dFitCh}(nn,:)]= getFrameCircleIntensity( green_frame, ...
      DcircleIndex{i},DringIndex{i}, subtract_bg);
    %red emission
    [data(i).intensity(nn,aFitCh) redTotalBG(i, nn) data(i).fitParam{aFitCh}(nn,:)]= getFrameCircleIntensity( red_frame, ...
      AcircleIndex{i},AringIndex{i}, subtract_bg);

  end
  
  if nn==1 
    txtOut = sprintf('%d\t/\t%d',nn,lastFrame);
    fprintf_r(txtOut,[],'reset');
  elseif nn==lastFrame || mod(nn,LOOPVAL)==0 
    txtOut = sprintf('%d\t/\t%d',nn,lastFrame);
    fprintf_r(txtOut,[]);
  end
  
end

data= subtractRingTimeBG( data, ringTime,  greenTotalBG, redTotalBG, firstGreenFrame,dFitCh,aFitCh,fretDataName);

fprintf('\nDone\n');

  %%%%%%%NESTED FUNCTION%%%%%%%%%%
  function  [circle_intensity total_background ringFitParams ]= getFrameCircleIntensity( frame , circleIndex, ringIndex, subtract_bg)
  %    [circle_intensity total_background ]= getFrameCircleIntensity( frame , circleIndex, ringIndex, subtract_bg)
  % outputs a vector circle_intensity which is the intensity for each image in the 
  % stack, for the given radius, and centre., with bg approximated as the intensity
  % of the ring which is subtacted from each point,*** IF bg subtraction is set***
  %    circle_intensity : vector size n, where n is the size of the 
  %          input stack
  %    stack    : input image
  %    centre : (x, y) vector, x, y, specifying the position of the
  %        point
  %    radius : integer number of pixels
  %    subtract_bg (true/false)

  %number of pixels in ring and circle
  totalRingPixels = numel(ringIndex);
  totalCirclePixels = numel(circleIndex);

  % calculate total intensities
  circle_intensity = sum(frame(circleIndex));
  circleAvg = circle_intensity/totalCirclePixels;
  ringAvg = sum(frame(ringIndex))./totalRingPixels;

  total_background = totalCirclePixels*ringAvg;
  if subtract_bg == true
    circle_intensity = circle_intensity - total_background;
  end
  
  %save all the per frame params
  ringFitParams = [circleAvg; ringAvg; totalCirclePixels; totalRingPixels];

  end
end

%-----------------------------------------------------------------------------------
function [ringIndex circleIndex] = calculateApertureIndex(frameSize, centre, innerCircleRadius, outerCircleRadius)
%   function [ringIndex circleIndex] = calculateApertureIndex(frameSize, centre, radius)
%
% calculate the points in index notation into the image frame which are inside the ring and
% the circle of given aperture size


if rem(innerCircleRadius , 1) ~= 0 
  error('the innerCircleRadius must be integer');
end

if rem(outerCircleRadius , 1) ~= 0 
  error('the outerCircleRadius must be integer');
end

x0 = centre(1);
y0 = centre(2);

% quantised circle - because we are using pixel coordinates
% the max makes sure the function only gives values greater than 1
circle = @(y,r, x0, y0) round( sqrt( r^2 - (y - y0).^2) + x0);

% boundaries of the outer and inner circle in x direction
xfinishO = circle( y0 - outerCircleRadius : y0 + outerCircleRadius, outerCircleRadius, x0, y0);
xstartO =  round(2*x0-xfinishO);
xfinishI = circle( y0 - innerCircleRadius : y0 + innerCircleRadius, innerCircleRadius, x0, y0);
xstartI =  round(2*x0-xfinishI);

% the (2*x0 - x) is just (x0-(x-x0)), ie the leftmost boundary of the circle
% (in case you were wondering)

ylim = frameSize(1);  % need to make sure we dont try to access out of 
xlim = frameSize(2);  % bounds values

% if j is out of bounds initialise it to an inbounds value
jstartO =  min( max( xstartO, 1) , xlim); 
jfinishO = min( max( xfinishO, 1), xlim);
jstartI =  min( max( xstartI, 1) , xlim); 
jfinishI = min( max( xfinishI, 1), xlim);

% sometimes jfinish is less than jstart due to rounding
% if this happens fix it
jfinishO(find(jfinishO<jstartO)) = jstartO(find(jfinishO<jstartO));
jfinishI(find(jfinishI<jstartI)) = jstartI(find(jfinishI<jstartI));

% boundaries in y direction
ystartO = round(y0)-outerCircleRadius; 
yfinishO = round(y0)+outerCircleRadius;
ystartI = round(y0)-innerCircleRadius; 
yfinishI = round(y0)+innerCircleRadius;

% if i is out of bounds initialise it to an inbounds value
istartO =  min( max( ystartO, 1), ylim); 
ifinishO = min( max( yfinishO, 1), ylim);
istartI =  min( max( ystartI, 1), ylim); 
ifinishI = min( max( yfinishI, 1), ylim);

% sometimes yfinish is less than ystart due to rounding
% if this happens fix it
ifinishO(find(ifinishO<istartO)) = istartO(find(ifinishO<istartO));
ifinishI(find(ifinishI<istartI)) = istartI(find(ifinishI<istartI));

% get the indexes of y-positions in the frame for the outer and inner circle
nO = 1;  
nI = 1;

% initialise the variables
maxSizeCircle = numel((round(y0)-outerCircleRadius : round(y0)+outerCircleRadius)) * numel(min(jstartO):max(jfinishO));
insideRing = zeros(maxSizeCircle , 2);
insideCircle = insideRing; 

%internal ring counter, m 
m = 0;
%internal circle counter, k
k = 0;

% loop over nO to create ring and circle
for i = istartO : ifinishO
  
  % left part of the ring
  if nO <= outerCircleRadius-innerCircleRadius
  
    RingPoints = numel(jstartO(nO):jfinishO(nO));
    
    if RingPoints > 0
      
      for s = 1 : RingPoints
      
      insideRing(m+s,1:2) =  [i, jstartO(nO)+s-1];
    
      end
      
      m = m + RingPoints;
      
    end
    
  end
  
  % middle part of the ring and the circle
  if nO > outerCircleRadius-innerCircleRadius
  if nO <= outerCircleRadius+innerCircleRadius+1
    
    leftRingPoints = numel(jstartO(nO):jstartI(nI));
    rightRingPoints = numel(jfinishI(nI):jfinishO(nO));
  
    if leftRingPoints + rightRingPoints > 0
      
      % create left part of the ring
      for s = 1 : leftRingPoints
      
        insideRing(m+s,1:2) =  [i, jstartO(nO)+s-1];
    
      end
      
      % create right part of the ring
      for s = 1 : rightRingPoints
    
	insideRing(m+leftRingPoints+s,1:2) =  [i, jfinishI(nI)+s-1];
  
      end
  
      m = m + leftRingPoints + rightRingPoints;
      
      % create circle inside the ring
      circlePoints = numel(jstartI(nI):jfinishI(nI))-2;
      
      if circlePoints > 0;
      
	for s = 1 : circlePoints
    
	  insideCircle(k+s,1:2) =  [i, jstartI(nI)+s];
  
	end
      
	k = k + circlePoints;
      
      end
    end
        
    %inner circle y-position index    
    nI = nI + 1;
  
  end
  end
  
  % right part of the ring
  if nO > outerCircleRadius + innerCircleRadius + 1
      
    RingPoints = numel(jstartO(nO):jfinishO(nO));
      
    if RingPoints > 0
          
      for s = 1 : RingPoints
          
        insideRing(m+s,1:2) =  [i, jstartO(nO)+s-1];
      
      end
          
      m = m + RingPoints;
          
    end
      
  end
  
  % outer circle y-position index
  nO = nO+1;

end

%convert the vectors from subscripts to indexes into frame
% 1: k-1 and m -1 because we dont want to try to index the zero padding values
ringIndex = sub2ind(frameSize, insideRing(1:m,1), insideRing(1:m,2));
circleIndex = sub2ind(frameSize, insideCircle(1:k,1), insideCircle(1:k,2));
end


%%-----------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------
function data= subtractRingTimeBG( data, ringTime,  greenTotalBG, redTotalBG, firstGreenFrame,dFitCh,aFitCh,fretDataName);
%function data= subtractRingTimeBG( data, ringTime,  greenTotalBG, redTotalBG, firstGreenFrame,dFitCh,aFitCh,fretDataName);
%  time averaged ring subtraction given a set of pre measured intensities and backgrounds
%
% ring time is the size of the forwards/backwards window ie half the total width.
%
% takes into account the alternation of the laser - splits into DexDem, AexAem, etc.

num_frames = size(greenTotalBG,2);
num_points = size(greenTotalBG,1);

if ringTime > num_frames/2
  error('Time average window is larger than total number of frames/ 2');
end

if firstGreenFrame == 0

  %pad the TotalBG matrices front and back for the sliding window;
  greenTotalBGPadded = [greenTotalBG(:,1:ringTime), greenTotalBG, greenTotalBG(:, end-ringTime+1:end)];
  redTotalBGPadded = [redTotalBG(:,1:ringTime), redTotalBG, redTotalBG(:, end-ringTime+1:end)];

  % initialise the TotalBGavg matrices
  greenTotalBGavg = zeros(size(greenTotalBG));
  redTotalBGavg = zeros(size(redTotalBG));

  %calculate the sliding window average
  for nn= 1:num_frames
    greenTotalBGavg(:,nn) = mean(greenTotalBGPadded(:,nn:nn+2*ringTime),2);
    redTotalBGavg(:,nn) = mean(redTotalBGPadded(:,nn:nn+2*ringTime),2);
  end
else
  %work out which one is the first red frame using a bit of simple modular arithmetic
  % nb matlab uses a modular system starting at zero so have to subtract 1 within the function 
  % and add one after the funtion.
  % its Mod(-2) because we also want the opposite frame
  firstRedFrame = mod(firstGreenFrame-2,2)+1;
  
  % split into the excitation channels
  DexDemTotalBG = greenTotalBG(:,firstGreenFrame:2:end);
  AexDemTotalBG = greenTotalBG(:,firstRedFrame:2:end);
  AexAemTotalBG = redTotalBG(:,firstRedFrame:2:end);
  DexAemTotalBG = redTotalBG(:,firstGreenFrame:2:end);
  
  numGreenFrames = size(DexDemTotalBG,2);
  numRedFrames = size(AexAemTotalBG,2);

  %pad the TotalBG matrices front and back for the sliding window;
  
  DexDemTotalBGpadded = [DexDemTotalBG(:,1:ringTime), DexDemTotalBG, DexDemTotalBG(:, end-ringTime+1:end)];
  AexDemTotalBGpadded = [AexDemTotalBG(:,1:ringTime), AexDemTotalBG, AexDemTotalBG(:, end-ringTime+1:end)];
  AexAemTotalBGpadded = [AexAemTotalBG(:,1:ringTime), AexAemTotalBG, AexAemTotalBG(:, end-ringTime+1:end)];
  DexAemTotalBGpadded = [DexAemTotalBG(:,1:ringTime), DexAemTotalBG, DexAemTotalBG(:, end-ringTime+1:end)];

  % initialise the TotalBGavg matrices

  DexDemTotalBGavg = zeros(size(DexDemTotalBG));
  AexDemTotalBGavg = zeros(size(AexDemTotalBG));
  AexAemTotalBGavg = zeros(size(AexAemTotalBG));
  DexAemTotalBGavg = zeros(size(DexAemTotalBG));

  greenTotalBGavg = zeros(size(greenTotalBG));
  redTotalBGavg = zeros(size(redTotalBG));
try  
  %calculate the sliding window average
  for nn = 1:numGreenFrames
    DexDemTotalBGavg(:,nn) = mean(DexDemTotalBGpadded(:,nn:nn+2*ringTime),2);
    DexAemTotalBGavg(:,nn) = mean(DexAemTotalBGpadded(:,nn:nn+2*ringTime),2);
  end

  for nn = 1:numRedFrames
    AexDemTotalBGavg(:,nn) = mean(AexDemTotalBGpadded(:,nn:nn+2*ringTime),2);
    AexAemTotalBGavg(:,nn) = mean(AexAemTotalBGpadded(:,nn:nn+2*ringTime),2);
  end
  % recombine the excitation channels
  greenTotalBGavg(:,firstGreenFrame:2:end) = DexDemTotalBGavg(:,:);
  greenTotalBGavg(:,firstRedFrame:2:end) = AexDemTotalBGavg(:,:);
  
  redTotalBGavg(:,firstRedFrame:2:end) = AexAemTotalBGavg(:,:);
  redTotalBGavg(:,firstGreenFrame:2:end) = DexAemTotalBGavg(:,:);
catch ME
  rethrow(ME);
end
end

for i = 1:num_points

  for nn = 1:num_frames
    data(i).intensity(nn,dFitCh) = data(i).intensity(nn,dFitCh) - greenTotalBGavg(i,nn);
    data(i).intensity(nn,aFitCh) = data(i).intensity(nn,aFitCh) - redTotalBGavg(i,nn);
    D = data(i).intensity(nn,dFitCh);
    A = data(i).intensity(nn,aFitCh);
    data(i).fretData(nn,:) = calculateE_CW(nn,D,A,fretDataName);
  end
end

end

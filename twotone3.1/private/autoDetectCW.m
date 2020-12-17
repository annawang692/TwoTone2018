function [twotoneDataOut, clusteredData, logText] = autoDetectCW(fileName,twotoneData);
% function [twotoneDataOut, clusteredData, logText] = autoDetectCW(fileName,twotoneData);
% perform autodetection on movies
% 
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 Alpha, released 101115
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%

% autoDetect the firstGreenFrame
imageLim = twotoneData.settings.imageSettings.channelImageLim;
alternationPeriod =1;
firstGreenFrame = 0;

tirfIm = TirfImage(fileName,firstGreenFrame, imageLim, alternationPeriod);
greenIm = getGreenStack(tirfIm);
numFrames = getNumFrames(greenIm);
imType = getImType(greenIm);
% calculate the averages 
% check the given limits are in bounds
avgFirst = twotoneData.settings.autoDetectSettings.averageFrameLim(1);
avgLast  = twotoneData.settings.autoDetectSettings.averageFrameLim(2);

if avgFirst > numFrames
  warning('twotone-autoDetectCmd:outOfBoundsFrame','Start frame to average is out of bounds, resetting to first frame');
  avgFirst = 1;
end

if avgLast > numFrames
  warning('twotone-autoDetectCmd:outOfBoundsFrame','Last frame to average is out of bounds, resetting to final frame');
  avgLast = numFrames;
end

[D,A] = calculateAveragesCW(tirfIm, avgFirst,avgLast);

% autoDetect points
Dch = find(strcmp(twotoneData.settings.imageSettings.aDetChannelName,'D'));
Ach = find(strcmp(twotoneData.settings.imageSettings.aDetChannelName,'A'));
thresholdD	= twotoneData.settings.autoDetectSettings.thresholds(Dch);
thresholdA	= twotoneData.settings.autoDetectSettings.thresholds(Ach);
bpDiscDiametre	= twotoneData.settings.autoDetectSettings.bandPassKernelDiameter;
windowSize	= twotoneData.settings.autoDetectSettings.fitSubImageRadius;

[pointsD pointsA] = ...
  autoDetectPoints(D,A,  ...
    thresholdD, thresholdA, bpDiscDiametre, windowSize);

paramNames = {'amplitude','sx','sy','eccentricity'};  
D = pointsD(:,1:2);
Dparams = pointsD(:,3:end);
A = pointsA(:,1:2);
Aparams = pointsA(:,3:end);

%check whether or not we are using a TFORM
applyTform  = twotoneData.settings.transformMatrixSettings.applyTform;
tformPath   = twotoneData.settings.transformMatrixSettings.fileName;
if applyTform==true
  try
    tformData = load(tformPath{1}); %Load the TFORM - loads it first as a structure
  catch ME %make the wrong tform crash sensible
    if strcmp(ME.identifier,'MATLAB:load:couldNotReadFile')
      error('autoDetectMain:TFORMerror','Unable to read TFORM file : No such file or directory.');
    else
      rethrow(ME);
    end
  end
  TFORM = tformData.TFORM;
else%do not use a tform 
  TFORM = [];
  tformData = [];
end

% link points
xmin = twotoneData.settings.autoDetectSettings.clusterDistanceThresh ;
filterChoice =twotoneData.settings.autoDetectSettings.linkageFilter;
[filteredClusters, allValidClusters,distanceDistribution]= associateCWchannels(D,A, xmin,Dparams,Aparams,paramNames,filterChoice,TFORM);

%filter the data

%get the filter parameter inputs
filterInputParams.applyNN	      = twotoneData.settings.autoDetectSettings.nearestNeighbor.apply;
filterInputParams.applyEccentricity   = twotoneData.settings.autoDetectSettings.ellipticity.apply;
filterInputParams.applySigma	      = twotoneData.settings.autoDetectSettings.PSFwidth.apply;
filterInputParams.NNlim		      = twotoneData.settings.autoDetectSettings.nearestNeighbor.thresh;
filterInputParams.eccLim	      = twotoneData.settings.autoDetectSettings.ellipticity.thresh;
filterInputParams.sigmaLim	      = twotoneData.settings.autoDetectSettings.PSFwidth.lim;

[filteredClusters] = filterLinkedDataCW(filteredClusters, filterInputParams);

[includedPos excludedPos]= getIncludedExcludedPos(filteredClusters,D,A,Dch,Ach);

% get the name of the movie
[path name extension] = fileparts(fileName);

%save the data
twotoneData.results.movieInfo.path  = path;
twotoneData.results.movieInfo.name  = fileName;
twotoneData.results.analysisDate    = datestr(now,'yymmd');
twotoneData.results.aDetImages{Dch}  = cast(D,imType);
twotoneData.results.aDetImages{Ach}  = cast(A,imType);
twotoneData.results.movieInfo.firstGreenFrame = firstGreenFrame;%save the alternation frame
twotoneData.results.aDetParamName = filteredClusters(1).paramNames;
%save the tform
twotoneData.results.transformMatrix.tform{1} = tformData;
%save all the position data
twotoneData.results.aDetMolPositions.included = includedPos;
twotoneData.results.aDetMolPositions.excluded = excludedPos;
%setup the final output s
%OUTPUT LOG FILE DATA
%get a bunch of log file parameters

% load the first green frame from the data
% get the number of detected particles in each channel
numD = size(D,1);
numA = size(A,1);
%get the number of linked particles
numLinkedParticles = numel(filteredClusters);

%print to screen
fprintf('First green frame %d\n',firstGreenFrame);
fprintf('Using TFORM: %s\n',tformPath{1});
fprintf('Particles in each channel:\n');
fprintf('D %d, A %d\n', numD,numA);
fprintf('Total linked particles: %d\n', numLinkedParticles);

%print to file
logText = '';
logText= sprintfappend(logText, 'First green frame %d\n',firstGreenFrame);
logText= sprintfappend(logText, 'Using TFORM: %s\n',tformPath{1});
logText= sprintfappend(logText, 'Particles in each channel:\n');
logText= sprintfappend(logText, 'D %d, A %d\n', numD,numA);
logText= sprintfappend(logText, 'Total linked particles: %d\n', numLinkedParticles);

%update the output
twotoneDataOut = twotoneData;
clusteredData = filteredClusters;
%------------------------------------------------------------------------------------------------
function [pointsD pointsA pointsDexAem] = ...
  autoDetectPoints(D,A, thresholdD, thresholdA, BPdiscDiametre,windowSize);
% function [pointsD pointsA ] = ...
%  autoDetectPoints(D,A,  thresholdD, thresholdA, BPdiscDiametre,windowSize);

% apply bandpass filter
[DFiltered,AFiltered ] = applyBPfilter(D,A, BPdiscDiametre);

%D
if isnan(thresholdD)
  pointsD = zeros(0,6);
else
  pointsD = findpeaks3(DFiltered,D, thresholdD, BPdiscDiametre,windowSize);
end
%A
if isnan(thresholdA)
  pointsA =zeros(0,6);
else
  pointsA = findpeaks3(AFiltered,A, thresholdA, BPdiscDiametre,windowSize);
end
%----------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------
function [DFiltered,AFiltered] = applyBPfilter(D,A, BPdiscDiametre)

%Check the disc diametre is valid
if isnan(BPdiscDiametre)... % is not a number
   || rem(BPdiscDiametre, 1)~=0 ... % is not integer
   || (BPdiscDiametre < 0)
  %dont update and thow an error
  success = 0;
  error('BPdisc must be integer >= 0');
else
  if BPdiscDiametre > 1
    DFiltered = bpass(D,1,BPdiscDiametre);
    AFiltered = bpass(A,1,BPdiscDiametre);
  else  % no low pass filtering only high pass for noise
    DFiltered = bpass(D,1,0);
    AFiltered = bpass(A,1,0);
  end

end
%----------------------------------------------------------------------------
function  [includedPos excludedPos]= getIncludedExcludedPos(filteredClusters,D,A, Dch,Ach);
%function [includedPos excludedPos]= getIncludedExcludedPos(filteredClusters,D,A, Dch,Ach);

% for the excludedPos include all the points initially, and strip the ones finally kept
excludedPos{Dch} = D;
excludedPos{Ach} = A;

nIncludedPos = numel(filteredClusters);
for i = 1:nIncludedPos
  %ad the particle to the included list
  includedPos{Dch}(i,:) = filteredClusters(i).d;
  %strip it from the excluded list
  xMatch = (excludedPos{Dch}(:,1)== filteredClusters(i).d(:,1));
  yMatch = (excludedPos{Dch}(:,2)== filteredClusters(i).d(:,2));
  matchRow = find(xMatch&yMatch);
  excludedPos{Dch}(matchRow,:)=[];

  %same for the other channels
  includedPos{Ach}(i,:) = filteredClusters(i).a;
  xMatch = (excludedPos{Ach}(:,1)== filteredClusters(i).a(:,1));
  yMatch = (excludedPos{Ach}(:,2)== filteredClusters(i).a(:,2));
  matchRow = find(xMatch&yMatch);
  excludedPos{Ach}(matchRow,:)=[];
end




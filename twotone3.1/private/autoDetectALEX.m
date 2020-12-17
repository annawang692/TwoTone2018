function [twotoneDataOut, clusteredData, logText] = autoDetectALEX(fileName,twotoneData);
% function [twotoneDataOut, clusteredData, logText] = autoDetectALEX(fileName,twotoneData);
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
alternationPeriod =twotoneData.settings.imageSettings.alternationPeriod;
firstGreenFrame = autoDetectALEXframe(fileName, imageLim);

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

[DexDem,AexAem, DexAem] = calculateAveragesALEX(tirfIm, firstGreenFrame, avgFirst,avgLast);

% autoDetect points
DDch = find(strcmp(twotoneData.settings.imageSettings.aDetChannelName,'DexDem'));
AAch = find(strcmp(twotoneData.settings.imageSettings.aDetChannelName,'AexAem'));
DAch = find(strcmp(twotoneData.settings.imageSettings.aDetChannelName,'DexAem'));
thresholdDD = twotoneData.settings.autoDetectSettings.thresholds(DDch);
thresholdAA = twotoneData.settings.autoDetectSettings.thresholds(AAch);
thresholdDA = twotoneData.settings.autoDetectSettings.thresholds(DAch);
bpDiscDiametre = twotoneData.settings.autoDetectSettings.bandPassKernelDiameter;
windowSize = twotoneData.settings.autoDetectSettings.fitSubImageRadius;

[pointsDexDem pointsAexAem pointsDexAem] = ...
  autoDetectPoints(DexDem,AexAem, DexAem, ...
    thresholdDD, thresholdAA, thresholdDA,  bpDiscDiametre,windowSize);

paramNames = {'amplitude','sx','sy','eccentricity'};  
DD = pointsDexDem(:,1:2);
DDparams = pointsDexDem(:,3:end);
AA = pointsAexAem(:,1:2);
AAparams = pointsAexAem(:,3:end);
DA = pointsDexAem(:,1:2);
DAparams = pointsDexAem(:,3:end);

%check whether or not we are using a TFORM
applyTform = twotoneData.settings.transformMatrixSettings.applyTform;
tformPath = twotoneData.settings.transformMatrixSettings.fileName;
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
[filteredClusters, allValidClusters,distanceDistribution]= associateALEXchannels(DD,AA,DA, xmin,DDparams,AAparams,DAparams,paramNames,filterChoice,TFORM);

%filter the data

%get the filter parameter inputs
filterInputParams.applyNN	      = twotoneData.settings.autoDetectSettings.nearestNeighbor.apply;
filterInputParams.applyEccentricity   = twotoneData.settings.autoDetectSettings.ellipticity.apply;
filterInputParams.applySigma	      = twotoneData.settings.autoDetectSettings.PSFwidth.apply;
filterInputParams.NNlim		      = twotoneData.settings.autoDetectSettings.nearestNeighbor.thresh;
filterInputParams.eccLim	      = twotoneData.settings.autoDetectSettings.ellipticity.thresh;
filterInputParams.sigmaLim	      = twotoneData.settings.autoDetectSettings.PSFwidth.lim;

[filteredClusters] = filterLinkedDataALEX(filteredClusters, filterInputParams);

[includedPos excludedPos]= getIncludedExcludedPos(filteredClusters,DD,AA,DA,DDch,AAch,DAch);

% get the name of the movie
[path name extension] = fileparts(fileName);

%save the data
twotoneData.results.movieInfo.path  = path;
twotoneData.results.movieInfo.name  = fileName;
twotoneData.results.analysisDate    = datestr(now,'yymmdd');
twotoneData.results.aDetImages{DDch}  = cast(DexDem,imType);
twotoneData.results.aDetImages{AAch}  = cast(AexAem,imType);
twotoneData.results.aDetImages{DAch}  = cast(DexAem,imType);
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
numDD = size(DD,1);
numAA = size(AA,1);
numDA = size(DA,1);
%get the number of linked particles
numLinkedParticles = numel(filteredClusters);

%print to screen
fprintf('First green frame %d\n',firstGreenFrame);
fprintf('Using TFORM: %s\n',tformPath{1});
fprintf('Particles in each channel:\n');
fprintf('DD %d, AA %d, DA %d\n', numDD,numAA,numDA);
fprintf('Total linked particles: %d\n', numLinkedParticles);

%print to file
logText = '';
logText= sprintfappend(logText, 'First green frame %d\n',firstGreenFrame);
logText= sprintfappend(logText, 'Using TFORM: %s\n',tformPath{1});
logText= sprintfappend(logText, 'Particles in each channel:\n');
logText= sprintfappend(logText, 'DD %d, AA %d, DA %d\n', numDD,numAA,numDA);
logText= sprintfappend(logText, 'Total linked particles: %d\n', numLinkedParticles);

%update the output
twotoneDataOut = twotoneData;
clusteredData = filteredClusters;
%------------------------------------------------------------------------------------------------
function [pointsDexDem pointsAexAem pointsDexAem] = ...
  autoDetectPoints(DexDem,AexAem, DexAem, thresholdDD, thresholdAA, thresholdDA, BPdiscDiametre,windowSize);
% function [pointsDexDem pointsAexAem pointsDexAem] = ...
%  autoDetectPoints(DexDem,AexAem, DexAem, thresholdDD, thresholdAA, thresholdDA, BPdiscDiametre,windowSize);

% apply bandpass filter
[DexDemFiltered,AexAemFiltered, DexAemFiltered] = applyBPfilter(DexDem,AexAem, DexAem,BPdiscDiametre);

%DD
if isnan(thresholdDD)
  pointsDexDem = zeros(0,6);
else
  pointsDexDem = findpeaks3(DexDemFiltered,DexDem, thresholdDD, BPdiscDiametre,windowSize);
end
%AA
if isnan(thresholdAA)
  pointsAexAem =zeros(0,6);
else
  pointsAexAem = findpeaks3(AexAemFiltered,AexAem, thresholdAA, BPdiscDiametre,windowSize);
end

%DA
if isnan(thresholdDA)
  pointsDexAem = zeros(0,6);
else
  pointsDexAem = findpeaks3(DexAemFiltered,DexAem, thresholdDA, BPdiscDiametre,windowSize);
end
%----------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------
function [DexDemFiltered,AexAemFiltered, DexAemFiltered] = applyBPfilter(DexDem,AexAem, DexAem,BPdiscDiametre)

%Check the disc diametre is valid
if isnan(BPdiscDiametre)... % is not a number
   || rem(BPdiscDiametre, 1)~=0 ... % is not integer
   || (BPdiscDiametre < 0)
  %dont update and thow an error
  success = 0;
  error('BPdisc must be integer >= 0');
else
  if BPdiscDiametre > 1
    DexDemFiltered = bpass(DexDem,1,BPdiscDiametre);
    AexAemFiltered = bpass(AexAem,1,BPdiscDiametre);
    DexAemFiltered = bpass(DexAem,1,BPdiscDiametre);
  else  % no low pass filtering only high pass for noise
    DexDemFiltered = bpass(DexDem,1,0);
    AexAemFiltered = bpass(AexAem,1,0);
    DexAemFiltered = bpass(DexAem,1,0);
  end

end
%----------------------------------------------------------------------------
function  [includedPos excludedPos]= getIncludedExcludedPos(filteredClusters,DD,AA,DA, DDch,AAch,DAch);
%function [includedPos excludedPos]= getIncludedExcludedPos(filteredClusters,DD,AA,DA, DDch,AAch,DAch);

% for the excludedPos include all the points initially, and strip the ones finally kept
excludedPos{DDch} = DD;
excludedPos{AAch} = AA;
excludedPos{DAch} = DA;

nIncludedPos = numel(filteredClusters);
for i = 1:nIncludedPos
  %add the particle to the included list
  includedPos{DDch}(i,:) = filteredClusters(i).dd;
  %strip it from the excluded list
  xMatch = (excludedPos{DDch}(:,1)== filteredClusters(i).dd(:,1));
  yMatch = (excludedPos{DDch}(:,2)== filteredClusters(i).dd(:,2));
  matchRow = find(xMatch&yMatch);
  excludedPos{DDch}(matchRow,:)=[];

  %same for the other channels
  includedPos{AAch}(i,:) = filteredClusters(i).aa;
  xMatch = (excludedPos{AAch}(:,1)== filteredClusters(i).aa(:,1));
  yMatch = (excludedPos{AAch}(:,2)== filteredClusters(i).aa(:,2));
  matchRow = find(xMatch&yMatch);
  excludedPos{AAch}(matchRow,:)=[];

  includedPos{DAch}(i,:) = filteredClusters(i).da;
  xMatch = (excludedPos{DAch}(:,1)== filteredClusters(i).da(:,1));
  yMatch = (excludedPos{DAch}(:,2)== filteredClusters(i).da(:,2));
  matchRow = find(xMatch&yMatch);
  excludedPos{DAch}(matchRow,:)=[];
end


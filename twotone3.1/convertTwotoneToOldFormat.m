function convertTwotoneToOldFormat(infilename,outfilename)
% convertTwotoneToOldFormat(infilename,outfilename)
%
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%
% FUNCTION: convertDataToOldTwotoneFormat
% DESCRIPTION:
%   Convert data from TwoTone versions >=3.1 to old format.
% INPUT: 
%   infilename  - Name of input data file
%   outfilename - Name to save output data file.
% OUTPUT: None

load(infilename);
[twotoneOldFormat] = convertDataToOldTwotoneFormat(twotoneData);
save(outfilename,'-struct','twotoneOldFormat');

%------------------------------------------------------------------------------
function  [twotoneOldFormat] = convertDataToOldTwotoneFormat(twotoneNewFormat)

twotoneOldFormat.analysis_info.path = twotoneNewFormat.results.movieInfo.name;%: original movie name, eg 'b16b17-g2r1-100ms_001.fits'
twotoneOldFormat.analysis_info.method = twotoneNewFormat.settings.fitSettings.algorithm;%: data extraction algorithm, eg 'fixedEllipseFit'
twotoneOldFormat.analysis_info.global_bg_subtract = false;%: IGNORE THIS, not used, I will delete it

    %*** analysis_info paramters from gaussian fit:
twotoneOldFormat.analysis_info.approx_radius = NaN;%: IGNORE THIS, not used, I will delete it
twotoneOldFormat.analysis_info.firstGreenFrame = twotoneNewFormat.results.movieInfo.firstGreenFrame; %: alternation (1/2), eg, 1

algorithm = twotoneNewFormat.settings.fitSettings.algorithm;
if strcmp(algorithm,'fixedGaussEllipse')|| strcmp(algorithm,'fixedGaussEllipse')|| strcmp(algorithm,'fixedGauss') || strcmp(algorithm,'freeGaussEllipse')
  %work out which row contains the WindowRadius
  windowRadiusRow = find(strcmp( twotoneNewFormat.settings.fitSettings.fitParamName, 'windowRadius'));
  fitMaxWidthRow  = find(strcmp( twotoneNewFormat.settings.fitSettings.fitParamName, 'maxWidth'));
  fitMinWidthRow  = find(strcmp( twotoneNewFormat.settings.fitSettings.fitParamName, 'minWidth'));
  twotoneOldFormat.analysis_info.FitCrop_radius = twotoneNewFormat.settings.fitSettings.param(windowRadiusRow) ;%: square window effective radius, integer, eg 6
  twotoneOldFormat.analysis_info.FitMaxwidth = twotoneNewFormat.settings.fitSettings.param(fitMaxWidthRow); %: 3
  twotoneOldFormat.analysis_info.FitMinwidth = twotoneNewFormat.settings.fitSettings.param(fitMinWidthRow);%: 1
elseif strcmp(algorithm,'ring')
  % *** analysis_info paramters from ring fit:
  insideRadiusRow = find(strcmp( twotoneNewFormat.settings.fitSettings.fitParamName, 'innerCircleRadius'));
  outsideRadiusRow = find(strcmp( twotoneNewFormat.settings.fitSettings.fitParamName, 'outerCircleRadius'));
  ringTimeRow = find(strcmp( twotoneNewFormat.settings.fitSettings.fitParamName, 'ringTime'));

  twotoneOldFormat.analysis_info.insideRadius = twotoneNewFormat.settings.fitSettings.param(insideRadiusRow);% : annulus inner radius
  twotoneOldFormat.analysis_info.outsideRadius = twotoneNewFormat.settings.fitSettings.param(outsideRadiusRow);% : annulus outer radius 
  twotoneOldFormat.analysis_info.intensityCircleRadius = twotoneNewFormat.settings.fitSettings.param(insideRadiusRow);% : circle (used to calculate raw intensity) radius
  twotoneOldFormat.analysis_info.intensityRingTime = twotoneNewFormat.settings.fitSettings.param(ringTimeRow);% : No of frames to include in calculation of BG from
  %annulus
else
  error('Unrecognised fitting algorithm');
end

nMolecules = numel(twotoneNewFormat.results.data);
if twotoneNewFormat.results.movieInfo.firstGreenFrame > 0
  dFitCh = find(strcmp(twotoneNewFormat.settings.imageSettings.channelName ,'Dem'));
  aFitCh = find(strcmp(twotoneNewFormat.settings.imageSettings.channelName ,'Aem'));
  ddAdetCh = find(strcmp(twotoneNewFormat.settings.imageSettings.aDetChannelName ,'DexDem'));
  aaAdetCh = find(strcmp(twotoneNewFormat.settings.imageSettings.aDetChannelName ,'AexAem'));
  daAdetCh = find(strcmp(twotoneNewFormat.settings.imageSettings.aDetChannelName ,'DexAem'));
else
  dFitCh = find(strcmp(twotoneNewFormat.settings.imageSettings.channelName ,'D'));
  aFitCh = find(strcmp(twotoneNewFormat.settings.imageSettings.channelName ,'A'));
  ddAdetCh = find(strcmp(twotoneNewFormat.settings.imageSettings.aDetChannelName ,'D'));
  aaAdetCh = find(strcmp(twotoneNewFormat.settings.imageSettings.aDetChannelName ,'A'));
  daAdetCh = find(strcmp(twotoneNewFormat.settings.imageSettings.aDetChannelName ,'D'));
end


for i = 1:nMolecules
  %Structure, 1 element for each analysed molecule
  twotoneOldFormat.intensities(i).DDadetPos = twotoneNewFormat.results.data(i).aDetData.aDetPos(ddAdetCh,:); %: [148.6474 44.3309] - autodetected position, DD
  twotoneOldFormat.intensities(i).AAadetPos = twotoneNewFormat.results.data(i).aDetData.aDetPos(aaAdetCh,:);%: [151.5059 44.3645] - autodetected position, AA - this is done for CW files too! (but its same as AA)
  twotoneOldFormat.intensities(i).DAadetPos = twotoneNewFormat.results.data(i).aDetData.aDetPos(daAdetCh,:);%: [151.4753 44.3734] - autodetected position, DA 
  twotoneOldFormat.intensities(i).DDadetParam = twotoneNewFormat.results.data(i).aDetData.aDetParam(ddAdetCh,:);%: [690.8210 1.3942 1.4818 0.3387] - autodetect parameters, defined in 'aDetParamNames' below, DD
  twotoneOldFormat.intensities(i).AAadetParam = twotoneNewFormat.results.data(i).aDetData.aDetParam(aaAdetCh,:);%: [1.0238e+03 1.5096 1.6669 0.4241]
  twotoneOldFormat.intensities(i).DAadetParam = twotoneNewFormat.results.data(i).aDetData.aDetParam(daAdetCh,:);%: [461.4777 1.4261 1.5661 0.4133]
  twotoneOldFormat.intensities(i).NNcur = twotoneNewFormat.results.data(i).aDetData.nearestNeighborDist;%: 8.8213 - nearest neighbour in all channels
  twotoneOldFormat.intensities(i).nDD = twotoneNewFormat.results.data(i).aDetData.nParticle(ddAdetCh);%: 1 -number of DD particles in cluster (0 /1 for any cluster
  twotoneOldFormat.intensities(i).nAA = twotoneNewFormat.results.data(i).aDetData.nParticle(ddAdetCh);%: 1 - as for DD
  twotoneOldFormat.intensities(i).nDA = twotoneNewFormat.results.data(i).aDetData.nParticle(ddAdetCh);%: 1 
  twotoneOldFormat.intensities(i).aDetParamNames = twotoneNewFormat.results.aDetParamName;%: {'amplitude'  'sx'  'sy'  'eccentricity'} - definintions of
  %params for DDadetPos etc
  twotoneOldFormat.intensities(i).fitParamNames = twotoneNewFormat.results.fitParamName;%: {'A0'  'sx'  'sy'  'BG'  'X'  'Y'  'theta'} - definintions
  %of params for greenFitParams etc
  twotoneOldFormat.intensities(i).green= twotoneNewFormat.results.data(i).intensity(:,dFitCh);%: [50x1 double] - [nFrame x 1] matrix
  twotoneOldFormat.intensities(i).red= twotoneNewFormat.results.data(i).intensity(:,aFitCh);%: [50x1 double]
  twotoneOldFormat.intensities(i).greenFitParams= twotoneNewFormat.results.data(i).fitParam{dFitCh}; %: [50x7 double] - [nFrame x nFitParam] matrix. For free
  %gaussian fitting, the XY position can change between each frame and this change
  %is recorded HERE, not in DDadetPos.
  twotoneOldFormat.intensities(i).redFitParams= twotoneNewFormat.results.data(i).fitParam{aFitCh}; %: [50x7 double] - as for green
  twotoneOldFormat.intensities(i).method= twotoneNewFormat.settings.fitSettings.algorithm; %: 'fixedEllipseFit' - fit method

  if strcmp(algorithm,'fixedGaussEllipse')|| strcmp(algorithm,'fixedGaussEllipse')|| strcmp(algorithm,'fixedGauss') || strcmp(algorithm,'freeGaussEllipse')
    twotoneOldFormat.intensities(i).crop_radius= twotoneNewFormat.settings.fitSettings.param(windowRadiusRow); %: 6 - square window effective radius
    twotoneOldFormat.intensities(i).maxwidth= twotoneNewFormat.settings.fitSettings.param(fitMaxWidthRow);%: 3 - as for analysis_info. ring fit will include the ring fit
    %specific params from analysis_info.
    twotoneOldFormat.intensities(i).minwidth= twotoneNewFormat.settings.fitSettings.param(fitMinWidthRow);%: 1 - as for analysis_info

  elseif strcmp(algorithm,'ring')
    % *** analysis_info paramters from ring fit:
    twotoneOldFormat.intensities(i).insideRadius = twotoneNewFormat.settings.fitSettings.param(insideRadiusRow);
    twotoneOldFormat.intensities(i).outsideRadius= twotoneNewFormat.settings.fitSettings.param(outsideRadiusRow);
    twotoneOldFormat.intensities(i).nBgFrameAvg =  twotoneNewFormat.settings.fitSettings.param(ringTimeRow);
    twotoneOldFormat.intensities(i).greenRingIndex = NaN; % this is not recorded anymore
    twotoneOldFormat.intensities(i).greenCircleIndex = NaN;
    twotoneOldFormat.intensities(i).redRingIndex = NaN;
    twotoneOldFormat.intensities(i).redCircleIndex = NaN;
  else
    error('Unrecognised fitting algorithm');
  end

  twotoneOldFormat.positionData.filteredClusters(i).NNcur = twotoneNewFormat.results.data(i).aDetData.nearestNeighborDist;%: 8.8213 
  twotoneOldFormat.positionData.filteredClusters(i).ndd = twotoneNewFormat.results.data(i).aDetData.nParticle(ddAdetCh);%: 1
  twotoneOldFormat.positionData.filteredClusters(i).naa = twotoneNewFormat.results.data(i).aDetData.nParticle(ddAdetCh);%: 1
  twotoneOldFormat.positionData.filteredClusters(i).nda = twotoneNewFormat.results.data(i).aDetData.nParticle(ddAdetCh);%: 1
  twotoneOldFormat.positionData.filteredClusters(i).dd = twotoneNewFormat.results.data(i).aDetData.aDetPos(ddAdetCh,:);%: [148.6474 44.3309]
  twotoneOldFormat.positionData.filteredClusters(i).aa = twotoneNewFormat.results.data(i).aDetData.aDetPos(aaAdetCh,:);%: [151.5059 44.3645]
  twotoneOldFormat.positionData.filteredClusters(i).da = twotoneNewFormat.results.data(i).aDetData.aDetPos(daAdetCh,:);%: [151.4753 44.3734]
  twotoneOldFormat.positionData.filteredClusters(i).ddparams = twotoneNewFormat.results.data(i).aDetData.aDetParam(ddAdetCh,:);%: [690.8210 1.3942 1.4818 0.3387]
  twotoneOldFormat.positionData.filteredClusters(i).aaparams = twotoneNewFormat.results.data(i).aDetData.aDetParam(aaAdetCh,:);%: [1.0238e+03 1.5096 1.6669 0.4241]
  twotoneOldFormat.positionData.filteredClusters(i).daparams = twotoneNewFormat.results.data(i).aDetData.aDetParam(daAdetCh,:);%: [461.4777 1.4261 1.5661 0.4133]
  twotoneOldFormat.positionData.filteredClusters(i).paramNames = twotoneNewFormat.results.aDetParamName;%: {'amplitude'  'sx'  'sy'  'eccentricity'}

end

twotoneOldFormat.positionData.movieName= twotoneNewFormat.results.movieInfo.name;%: 'b16b17-g2r1-100ms_001.fits' - movie name
twotoneOldFormat.positionData.moviePath= fullfile(twotoneNewFormat.results.movieInfo.path,twotoneNewFormat.results.movieInfo.name);%: 'b16b17-g2r1-100ms_001.fits' - full path to movie if not in the local directory
twotoneOldFormat.positionData.filter = twotoneNewFormat.settings.autoDetectSettings.linkageFilter;%: {'DexDem&&AexAem'} - linkage filter applied between channels
twotoneOldFormat.positionData.DexDem = twotoneNewFormat.results.aDetImages{ddAdetCh};%: [258x256 uint16] - autodetected image
twotoneOldFormat.positionData.AexAem = twotoneNewFormat.results.aDetImages{aaAdetCh};%: [258x256 uint16] - autodetected image
twotoneOldFormat.positionData.DexAem = twotoneNewFormat.results.aDetImages{daAdetCh};%: [258x256 uint16] - autodetected image
twotoneOldFormat.positionData.fitsHeader = twotoneNewFormat.results.movieInfo.imageInfo;%: [1x1 struct] - header of the input fits file
twotoneOldFormat.positionData.firstGreenFrame = twotoneNewFormat.results.movieInfo.firstGreenFrame;%: 1 - determine alternation (1/2)
twotoneOldFormat.positionData.TFORM = twotoneNewFormat.results.transformMatrix.tform{1};%: [1x1 struct] - transformation matrix used to map red positions to coordinate sys of green channel.
twotoneOldFormat.positionData.TFORMfilename = twotoneNewFormat.settings.transformMatrixSettings.fileName;%: 'b16b17LEAK.tform.mat' - image used to generate transform
includedPos = twotoneNewFormat.results.aDetMolPositions.included;
excludedPos = twotoneNewFormat.results.aDetMolPositions.excluded;

twotoneOldFormat.positionData.DexDemPos = [includedPos{ddAdetCh}; excludedPos{ddAdetCh} ] ;%: [160x6 double] - all detected DD positions
twotoneOldFormat.positionData.DexAemPos = [includedPos{daAdetCh}; excludedPos{daAdetCh} ] ;%: [148x6 double] - all detected AA positions
twotoneOldFormat.positionData.AexAemPos = [includedPos{aaAdetCh}; excludedPos{aaAdetCh} ]; %: [152x6 double] - all detected DA positions
twotoneOldFormat.positionData.applyNN = twotoneNewFormat.settings.autoDetectSettings.nearestNeighbor.apply;%: 1 - use nearest neighbour threshold? (0/1)
twotoneOldFormat.positionData.applyEccentricity = twotoneNewFormat.settings.autoDetectSettings.ellipticity.apply;%: 0 - use eccentricity threshold? (0/1)
twotoneOldFormat.positionData.applySigma = twotoneNewFormat.settings.autoDetectSettings.PSFwidth.apply;%: 0 - use sigma threshold? (0/1)
twotoneOldFormat.positionData.NNlim = twotoneNewFormat.settings.autoDetectSettings.nearestNeighbor.thresh;%: 8.5000 - min nearest neighbour
twotoneOldFormat.positionData.eccLim = twotoneNewFormat.settings.autoDetectSettings.ellipticity.thresh;%: [0 Inf] - eccentricity min/max
twotoneOldFormat.positionData.sigmaLim = twotoneNewFormat.settings.autoDetectSettings.PSFwidth.lim;%: [0 Inf] - autodetction gaussian fit width min/max


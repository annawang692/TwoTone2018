function initializeTwotoneSetting(specifyFile);
% initialize twotone default settings;
%   specifyFile = 'ALEXonly', 'CWonly', 'algSettingOnly'(default = all);
% 
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 Alpha, released 101115
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%

if ~exist('specifyFile', 'var') || strcmp(specifyFile,'algSettingOnly')

  twotoneGeneralSetting.linkageFilterOpts_ALEX={'showAll', 'DexDem', 'AexAem','DexAem',...
      'DexDem&&AexAem', 'DexDem&&DexAem', 'AexAem&&DexAem', 'DexDem&&AexAem&&DexAem',...
      'DexDem||AexAem||DexAem', 'DexDem||AexAem'};
  twotoneGeneralSetting.linkageFilterOpts_CW={'showAll', 'D', 'A','D&&A','D||A'};
  % Fitting algorithm
  twotoneGeneralSetting.analysisAlgorithm={'fixedGaussEllipse','ring','fixedGauss','freeGaussEllipse'};
  %fixedGaussEllipse;
  % 2D elliptical fit, fixed position;
  twotoneGeneralSetting.algorithmParam{1}={'windowRadius','minWidth','maxWidth'} ;
  twotoneGeneralSetting.algorithmDefaultVal{1} = [6 1 2];
  %ring;
  %Ring fit (Aperture photometry);
  twotoneGeneralSetting.algorithmParam{2}={'innerCircleRadius','outerCircleRadius','ringTime'};
  twotoneGeneralSetting.algorithmDefaultVal{2} = [4 6 3];

  %fixedgauss;
  % 2D circular Gaussian fit, fixed position;
  twotoneGeneralSetting.algorithmParam{3}={'windowRadius','minWidth','maxWidth'} ;
  twotoneGeneralSetting.algorithmDefaultVal{3} = [6 1 2];

  %freegaussellipse;
  % 2D elliptical Gaussian fit, free position;
  twotoneGeneralSetting.algorithmParam{4}={'windowRadius','minWidth','maxWidth','maxPositionChange'} ;
  twotoneGeneralSetting.algorithmDefaultVal{4} = [6 1 2 3];
  
  twotoneGeneralSetting.PSFwidthEstimate = 1.5;
  
  twotoneGeneralSetting.outputOldTwotoneFormat = false;
  twotoneGeneralSetting.outputTwotoneTextData  = false;

  save(twotoneInstallDir('','twotoneGeneralSetting.mat'),  'twotoneGeneralSetting');
end

if ~exist('specifyFile', 'var') || strcmp(specifyFile,'ALEXonly')
  clear('twotoneData');
  %twotoneSettings_ALEX;

  twotoneData.settings.fileAppendString = '.fitResult.mat';

  twotoneData.settings.imageSettings.twotoneVersion = '3.1.0alpha';
  twotoneData.settings.imageSettings.alternationPeriod = 2;
  twotoneData.settings.imageSettings.nChannel = 2;
  twotoneData.settings.imageSettings.channelName = {'Dem','Aem'};
  twotoneData.settings.imageSettings.channelImageLim= [[1,256,1,256];[257,512,1,256]];
  twotoneData.settings.imageSettings.nADetChannel= 3;
  twotoneData.settings.imageSettings.aDetChannelName={'DexDem', 'AexAem', 'DexAem' };

  twotoneData.settings.transformMatrixSettings.applyTform = 0;
  twotoneData.settings.transformMatrixSettings.channelLinkage= {[1,2]}; %specifies which channel the transform is from and to (in this case 1->2) ;
  twotoneData.settings.transformMatrixSettings.fileName= {''};

  twotoneData.settings.fitSettings.algorithm= 'fixedGauss';
  twotoneData.settings.fitSettings.fitParamName= {'windowRadius','minWidth','maxWidth'};
  twotoneData.settings.fitSettings.param= [6, 1,2];

  twotoneData.settings.fretDataName = {'t_Dex','t_Aex','DD','DA','AD','AA','E','S'};
  

  twotoneData.settings.autoDetectSettings.averageFrameLim= [2, 10];
  twotoneData.settings.autoDetectSettings.bandPassKernelDiameter= 3;
  twotoneData.settings.autoDetectSettings.fitSubImageRadius= 3;
  twotoneData.settings.autoDetectSettings.thresholds =  [15,15, 10];
  twotoneData.settings.autoDetectSettings.clusterDistanceThresh  = 3;
  twotoneData.settings.autoDetectSettings.nearestNeighbor.apply  = 0;
  twotoneData.settings.autoDetectSettings.nearestNeighbor.thresh = 6;
  twotoneData.settings.autoDetectSettings.ellipticity.apply  = 0;
  twotoneData.settings.autoDetectSettings.ellipticity.thresh = 0.6;
  twotoneData.settings.autoDetectSettings.PSFwidth.apply= 0;
  twotoneData.settings.autoDetectSettings.PSFwidth.lim= [1, 2];
  twotoneData.settings.autoDetectSettings.linkageFilter  = 'DexDem&&AexAem';

  save(twotoneInstallDir('','twotoneDefaultSettings_ALEX.mat'),  'twotoneData');
end

if ~exist('specifyFile', 'var') || strcmp(specifyFile,'CWonly')
  %twotoneSettings_CW;
  clear('twotoneData');
  twotoneData.settings.fileAppendString = '.fitResult.mat';

  twotoneData.settings.imageSettings.twotoneVersion = '3.1.0alpha';
  twotoneData.settings.imageSettings.alternationPeriod = 1;
  twotoneData.settings.imageSettings.nChannel = 2;
  twotoneData.settings.imageSettings.channelName = {'D','A'};
  twotoneData.settings.imageSettings.channelImageLim= [[1,256,1,256];[257,512,1,256]];
  twotoneData.settings.imageSettings.nADetChannel= 2;
  twotoneData.settings.imageSettings.aDetChannelName={'D', 'A'};

  twotoneData.settings.transformMatrixSettings.applyTform = 0;
  twotoneData.settings.transformMatrixSettings.channelLinkage= {[1,2]}; %specifies which channel the transform is from and to (in this case 1->2) ;
  twotoneData.settings.transformMatrixSettings.fileName= {''};

  twotoneData.settings.fitSettings.algorithm= 'fixedGauss';
  twotoneData.settings.fitSettings.fitParamName= {'windowRadius','minWidth','maxWidth'};
  twotoneData.settings.fitSettings.param= [6, 1,2];

  twotoneData.settings.fretDataName = {'t','D','A','E'};

  twotoneData.settings.autoDetectSettings.averageFrameLim= [2, 10];
  twotoneData.settings.autoDetectSettings.bandPassKernelDiameter= 3;
  twotoneData.settings.autoDetectSettings.fitSubImageRadius= 3;
  twotoneData.settings.autoDetectSettings.thresholds=  [15,15];
  twotoneData.settings.autoDetectSettings.clusterDistanceThresh= 3;
  twotoneData.settings.autoDetectSettings.nearestNeighbor.apply= 0;
  twotoneData.settings.autoDetectSettings.nearestNeighbor.thresh= 6;
  twotoneData.settings.autoDetectSettings.ellipticity.apply= 0;
  twotoneData.settings.autoDetectSettings.ellipticity.thresh= 0.6;
  twotoneData.settings.autoDetectSettings.PSFwidth.apply= 0;
  twotoneData.settings.autoDetectSettings.PSFwidth.lim= [1, 2];
  twotoneData.settings.autoDetectSettings.linkageFilter= 'D';

  save(twotoneInstallDir('','twotoneDefaultSettings_CW.mat'), 'twotoneData');

end

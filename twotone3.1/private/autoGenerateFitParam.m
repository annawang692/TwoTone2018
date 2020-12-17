function [twotoneData, twotoneGeneralSetting]= autoGenerateFitParam(twotoneData, twotoneGeneralSetting,width)
% function autoParamOut = autoGenerateFitParam(width)
%  calculate default values for the paramters based on 
%  the size of the PSF
% 
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 Alpha, released 101115
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%
autoParam = estimateParam(width);
[twotoneData, twotoneGeneralSetting] = updateParam(twotoneData, twotoneGeneralSetting,autoParam );


%-----------------------------------------------
function autoParam = estimateParam(width)
% calculate the appropriate settings

%autoParam.autoDetWindowRadius = max(2, round(2*width));
%autoParam.profileFittingRadius = max(3, round(4*width));
%110419 changed to improve 2x2 binning performance
autoParam.autoDetWindowRadius = max(3, round(2*width));
autoParam.profileFittingRadius = max(4, round(4*width));

autoParam.profileFittingMaxPosChange = autoParam.profileFittingRadius/2;
fitMinWidth = max(0,width - 0.5);
fitMaxWidth = width + 0.5;
autoParam.widthLim = [fitMinWidth fitMaxWidth];
autoParam.nnLim = 4*width;

autoParam.aperturePhotInner = floor(3*width);
autoParam.aperturePhotOuter = round(sqrt(2)*3*width); %the sqrt 2 factor gives roughly equal areas for the inner and outer radii
if autoParam.aperturePhotInner == autoParam.aperturePhotOuter
  autoParam.aperturePhotOuter = autoParam.aperturePhotOuter + 1;
end

autoParam.aperturePhotFrameAvg = 3; %default setting

%autoParam.bpasskerneldiametre = max(round(2*width),1);
%110419 changed to improve 2x2 binning performance
autoParam.bpasskerneldiametre = max(round(2*width),3);

%-----------------------------------------------
function [twotoneData, twotoneGeneralSetting] = updateParam(twotoneData, twotoneGeneralSetting,autoParam)


%update twotoneGeneralSetting first

%update each algorithms parameter set
%update fixedGaussEllipse
currentAlg='fixedGaussEllipse';
cellNo = find(strcmp(twotoneGeneralSetting.analysisAlgorithm,currentAlg));
currentParam = 'windowRadius';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) = autoParam.profileFittingRadius;
currentParam = 'minWidth';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) = autoParam.widthLim(1);
currentParam = 'maxWidth';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) = autoParam.widthLim(2);

%update fixedGauss
currentAlg='fixedGauss';
cellNo = find(strcmp(twotoneGeneralSetting.analysisAlgorithm,currentAlg));
currentParam = 'windowRadius';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) = autoParam.profileFittingRadius;
currentParam = 'minWidth';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) = autoParam.widthLim(1);
currentParam = 'maxWidth';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) = autoParam.widthLim(2);

%update freeGaussEllipse
currentAlg='freeGaussEllipse';
cellNo = find(strcmp(twotoneGeneralSetting.analysisAlgorithm,currentAlg));
currentParam = 'windowRadius';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) = autoParam.profileFittingRadius;
currentParam = 'minWidth';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) = autoParam.widthLim(1);
currentParam = 'maxWidth';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) = autoParam.widthLim(2);
currentParam = 'maxPositionChange';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) = autoParam.profileFittingMaxPosChange;

%update fixedGauss
currentAlg='ring';
cellNo = find(strcmp(twotoneGeneralSetting.analysisAlgorithm,currentAlg));
currentParam = 'innerCircleRadius';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) =autoParam.aperturePhotInner;
currentParam = 'outerCircleRadius';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) = autoParam.aperturePhotOuter;
currentParam = 'ringTime';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) =autoParam.aperturePhotFrameAvg;

% now update twotoneData fit algorithm settings
selectedAlgNo = find(strcmp(twotoneGeneralSetting.analysisAlgorithm,twotoneData.settings.fitSettings.algorithm));
twotoneData.settings.fitSettings.param=twotoneGeneralSetting.algorithmDefaultVal{selectedAlgNo};

%now update the other twotoneData settings
twotoneData.settings.autoDetectSettings.fitSubImageRadius	= autoParam.autoDetWindowRadius ;
twotoneData.settings.autoDetectSettings.bandPassKernelDiameter	= autoParam.bpasskerneldiametre ;
twotoneData.settings.autoDetectSettings.nearestNeighbor.thresh	= autoParam.nnLim;
twotoneData.settings.autoDetectSettings.PSFwidth.lim		= autoParam.widthLim;

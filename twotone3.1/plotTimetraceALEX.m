function [absoluteIndex plotDataOut] = plotTimetraceALEX(plotData,varargin)
% [absoluteIndex plotDataOut]= plotTimetraceALEX(plotData, varargin)
% 
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 , released 110426
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%
% FUNCTION: plotTimetraceALEX
% DESCRIPTION:
%   Plot fret, stoichiometry and intensity timetraces from TwoTone output data. Perform filtering on the data.
% 
% INPUTS:
%   plotData - 3 options:
%	'file1.mat' - name of single file for analysis
%	{'file1.mat', 'file2.mat',...,'filen.mat'} - cell containing the names of mulitple files for analysis
%	plotDataMatrix - the matrix output from plotHistALEX can be re-input (eg to select a different particle) - in this case filtering arguments will be IGNORED!
%   filter (optional): (any number of these may be applied): {'filterType','filterParam','>' OR '<', Value}
%	FilterType:  'fretData', 
%	  FilterParam:   't_Dex', 't_Aex', 'DD', 'DA','AA', 'AA', 'E', 'S'
%	FilterType: 'all_aDetParam'
%	  FilterParam:  'nearestNeighbor'
%	FilterType:    'DD_aDetParam'
%	FilterType:    'DA_aDetParam'
%	FilterType:    'AA_aDetParam'
%	  FilterParam:     'amplitude', 'sx', 'sy', 'eccentricity', 
%	  
%	FilterType:    'DD_aDetPos'
%	FilterType:    'DA_aDetPos'
%	FilterType:    'AA_aDetPos'
%	  FilterParam:     'X' 'Y'
%	Range: '<' or '>'
%	Value: Scalar number
%   'RelativeIndex', relativeIndex (optional) - n-th particle which meets thresholds in dataset. Multiple particles may be specified
%   'AbsoluteIndex', [movieNo, ParticleNo] (optional) - n-th particle in m-th movie - multiple particles may be specified on separate rows.
% 'IntegrationTime', integrationTime (optional) - Specify the duration of a single frame of the raw data - used to plot the time axis. If not specified an integration time of 1 is used.
% 
% OUTPUTS: 
%   absoluteIndex: [movieNo, ParticleNo] - vector of all plotted particles
%   plotDataOut: structure containing plotted data.
%     Contains:
%	plotOutput.dataName ={'movieNo' 'particleNo' 'DDX' 'DDY' 't_Dex' 't_Aex' 'DD' 'DA' 'AD' 'AA' 'E' 'S'}
%	plotOutput.data - matrix containing columns specified in dataName. Each row is one frame

n = numel(varargin);
i = 1;
integrationTime = 1;
filterArg = {};
selectParticleArg = {};
while i <= n
  if strcmp(varargin{i},'IntegrationTime')
    integrationTime = varargin{i+1};
    i=i+2;
  elseif iscell(varargin{i})
    filterArg = {filterArg{:}, varargin{:}};
    i = i + 1;
  else
    selectParticleArg= {selectParticleArg{:},varargin{i}};
    i = i + 1;
  end
end

if isstr(plotData) || iscell(plotData)
  plotData = filterData_ALEX(plotData, filterArg{:});
end

if numel(selectParticleArg) == 0
  selectParticleArg = {'RelativeIndex',1};
end 

[plotData absoluteIndex] = selectParticle(plotData,selectParticleArg{:});

clf; fig1 = gcf; 
set(fig1,'color','w','name','FRET and Intensity Timetrace')

%time series Fret
ax(1) = subplot(2,1,1);
title(['Movie ',num2str(absoluteIndex(:,1)'),', Particle ', num2str(absoluteIndex(:,2)')]);
movieCol	  = find(strcmp(plotData.dataName,'movieNo'));
particleCol	  = find(strcmp(plotData.dataName,'particleNo'));
tDexCol	  = find(strcmp(plotData.dataName,'t_Dex'));
tAexCol	  = find(strcmp(plotData.dataName,'t_Aex'));
eCol	  = find(strcmp(plotData.dataName,'E'));
sCol	  = find(strcmp(plotData.dataName,'S'));
DDCol	  = find(strcmp(plotData.dataName,'DD'));
AACol	  = find(strcmp(plotData.dataName,'AA'));
DACol	  = find(strcmp(plotData.dataName,'DA'));

t_Dex =  plotData.data(:,tDexCol)*integrationTime;
t_Aex =  plotData.data(:,tAexCol)*integrationTime;
t_S = mean([t_Dex,t_Aex],2);
E = plotData.data(:,eCol);
S = plotData.data(:,sCol);
xlabel('t');
ylabel('\color{red} E\color{black}/ S');

hold on;
stairs(t_Dex,E, 'r');
stairs(t_S,S, 'k');
ylim([0 1]);
set(gca,'YGrid','on');

maxVal = max([E(:); S(:)]);
minVal = min([E(:); S(:)]);

%stairs Intensities
ax(2) = subplot(2,1,2);

DD = plotData.data(:,DDCol);
AA = plotData.data(:,AACol);
DA = plotData.data(:,DACol);

D = DD + DA;

hold on;
stairs(t_Dex,DD, 'g');
stairs(t_Dex,DA, 'r');
stairs(t_Dex,D , 'b');
stairs(t_Aex,AA, 'k');
set(gca,'YGrid','on');
ylabel('\color{green}DD/\color{red} DA/\color{blue} Dex/ \color{black} AA (au)','HorizontalAlignment','Center');
xlabel('t');

maxIntensity = max([DD(:); DA(:); D(:); AA(:)]);
minIntensity = min([DD(:); DA(:); D(:); AA(:)]);
ylim([min(minIntensity,0) 1.2*maxIntensity]);
%link the x-axes of the two subplots
linkaxes(ax,'x');

plotDataOut = plotData;

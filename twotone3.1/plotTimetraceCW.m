function [absoluteIndex plotDataOut]= plotTimetraceCW(plotData,varargin)
% [absoluteIndex plotDataOut]= plotTimetraceCW(plotData, varargin)
% 
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 , released 110426
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%
% FUNCTION: plotTimetraceCW
%
% DESCRIPTION:
%   Plot fret, and intensity timetraces from TwoTone output data. Perform filtering on the data.
% 
% INPUTS:
%   plotData - 3 options:
%	'file1.mat' - name of single file for analysis
%	{'file1.mat', 'file2.mat',...,'filen.mat'} - cell containing the names of mulitple files for analysis
%	plotDataMatrix - the matrix output from plotHistALEX can be re-input (eg to select a different particle) - in this case filtering arguments will be IGNORED!
%   filter (optional): (any number of these may be applied): {'filterType','filterParam','>' OR '<', Value}
%	FilterType:  'fretData', 
%	  FilterParam:   't', 'D','A','E'
%	FilterType: 'all_aDetParam'
%	  FilterParam:  'nearestNeighbor'
%	FilterType:    'D_aDetParam'
%	FilterType:    'A_aDetParam'
%	  FilterParam:     'amplitude', 'sx', 'sy', 'eccentricity', 
%	FilterType:    'D_aDetPos'
%	FilterType:    'A_aDetPos'
%	  FilterParam:     'X' 'Y'
%	Range: '<' or '>'
%	Value: Scalar number
%   'NBin', nBin (optional) - number of bins in histogram
%   'RelativeIndex', relativeIndex (optional) - n-th particle which meets thresholds in dataset. Multiple particles may be specified
%   'AbsoluteIndex', [movieNo, ParticleNo] (optional) - n-th particle in m-th movie - multiple particles may be specified on separate rows.
% 'IntegrationTime', integrationTime (optional) - Specify the duration of a single frame of the raw data - used to plot the time axis. If not specified an integration time of 1 is used.
% 
% OUTPUTS: 
%   absoluteIndex: [movieNo, ParticleNo] - vector of all plotted particles
%   plotDataOut: structure containing plotted data.
%     Contains:
%	plotOutput.dataName ={'movieNo' 'particleNo' 'DX' 'DY' 't' 'D' 'A' 'E'}
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
  plotData = filterData_CW(plotData, filterArg{:});
end

if numel(selectParticleArg) == 0
  selectParticleArg= {'RelativeIndex',1};
end

[plotData absoluteIndex] = selectParticle(plotData,selectParticleArg{:});

clf; fig1 = gcf; 
set(fig1,'color','w','name','FRET and Intensity Timetrace')

%time series Fret
ax(1) = subplot(2,1,1);
title(['Movie ',num2str(absoluteIndex(:,1)'),', Particle ', num2str(absoluteIndex(:,2)')]);
movieCol	  = find(strcmp(plotData.dataName,'movieNo'));
particleCol	  = find(strcmp(plotData.dataName,'particleNo'));
tCol	  = find(strcmp(plotData.dataName,'t'));
eCol	  = find(strcmp(plotData.dataName,'E'));
DCol	  = find(strcmp(plotData.dataName,'D'));
ACol	  = find(strcmp(plotData.dataName,'A'));

t =  plotData.data(:,tCol)*integrationTime;
E = plotData.data(:,eCol);
xlabel('t');
ylabel('\color{red} E');

hold on;
stairs(t,E, 'r');
ylim([0 1]);
set(gca,'YGrid','on');

maxVal = max(E);
minVal = min(E);

%stairs Intensities
ax(2) = subplot(2,1,2);

D = plotData.data(:,DCol);
A = plotData.data(:,ACol);
N = D + A;

hold on;
stairs(t,D, 'g');
stairs(t,A, 'r');
stairs(t,N , 'b');
set(gca,'YGrid','on');
ylabel('\color{green}D/\color{red} A/\color{blue} N (au)','HorizontalAlignment','Center');
xlabel('t');

maxIntensity = max([D(:); A(:); N(:)]);
minIntensity = min([D(:); A(:); N(:)]);
ylim([min(minIntensity,0) 1.2*maxIntensity]);
%link the x-axes of the two subplots
linkaxes(ax,'x');


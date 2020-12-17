function [absoluteIndex plotDataOut Ex En]= plotHistCW(plotData, varargin)
% [absoluteIndex plotDataOut Ex En]= plotHistCW(plotData, varargin)
% 
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0, released 110426
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%
% FUNCTION:  plotHistCW
% DESCRIPTION:
%   Plot E histogram from TwoTone output data. Perform filtering on the data.
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
%   ELim, [eMin eMax] (optional) - maximum and minimum fret values
%   'RelativeIndex', relativeIndex (optional) - n-th particle which meets thresholds in dataset. Multiple particles may be specified
%   'AbsoluteIndex', [movieNo, ParticleNo] (optional) - n-th particle in m-th movie - multiple particles may be specified on separate rows.
% 
% OUTPUTS: 
%   absoluteIndex: [movieNo, ParticleNo] - vector of all plotted particles
%   plotDataOut: structure containing plotted data.
%     Contains:
%	plotOutput.dataName ={'movieNo' 'particleNo' 'DX' 'DY' 't' 'D' 'A' 'E'}
%	plotOutput.data - matrix containing columns specified in dataName. Each row is one frame
%   Ex, En - fret histogram data

n = numel(varargin);
i = 1;
nBin= 100;
eLim = [0 1];
filterArg = {};
selectParticleArg = {};
while i <= n
  if strcmp(varargin{i},'NBin')
    nBin = varargin{i+1};
    i=i+2;
  elseif strcmp(varargin{i},'ELim')
    eLim = varargin{i+1};
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

[plotData absoluteIndex] = selectParticle(plotData,selectParticleArg{:});

eCol = find(strcmp(plotData.dataName,'E'));

Ematrix = plotData.data(:,eCol);

%calculate the limits for the plots
Min_E = eLim(1);
Max_E = eLim(2);

Ematrix(find(Ematrix(:,1)<Min_E | Ematrix(:,1)>Max_E)) = []; 

%make sure we get exactly 100 bins by using bin centre
[eCentre eEdge]= bincentre(Min_E,Max_E,nBin);

%plot it here

[TopN, TopX]= hist(Ematrix(:,1),eCentre);
hE = bar(TopX,TopN);
Ex = TopX;
En = TopN;
ax1=gca;
xlabel('E','FontSize',14,'Parent',ax1);
ylabel('Freq.','FontSize',14,'Parent',ax1);
h01 = findobj(gca,'Type','patch');
set(h01,'FaceColor','b','EdgeColor','b');

% set the axes labels and limits
axis(ax1,[Min_E Max_E 0 max(TopN)*1.2]);

eTickLabels = [Min_E:((Max_E - Min_E)/5):Max_E];

plotDataOut = plotData;

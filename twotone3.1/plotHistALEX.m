function [absoluteIndex plotDataOut Ex En Sx Sn]= plotHistALEX(plotData, varargin)
% [absoluteIndex plotDataOut Ex En Sx Sn]= plotHistALEX(plotData, varargin)
%
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 , released 110426
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%
% FUNCTION:  plotHistALEX
% DESCRIPTION:
%   Plot ES histogram from TwoTone output data. Perform filtering on the data.
% 
% INPUTS:
%   plotData - 3 options:
%	'file1.mat' - name of single file for analysis
%	{'file1.mat', 'file2.mat',...,'filen.mat'} - cell containing the names of mulitple files for analysis
%	plotDataMatrix - the matrix output from plotHistALEX can be re-input (eg to select a different particle)
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
%   'NBin', nBin (optional) - number of bins in histogram
%   ELim, [eMin eMax] (optional) - maximum and minimum fret values
%   SLim, [sMin sMax] (optional) - max and min stoichiometry values
%   'RelativeIndex', relativeIndex (optional) - n-th particle which meets thresholds in dataset. Multiple particles may be specified
%   'AbsoluteIndex', [movieNo, ParticleNo] (optional) - n-th particle in m-th movie - multiple particles may be specified on separate rows.
% 
% OUTPUTS: 
%   absoluteIndex: [movieNo, ParticleNo] - vector of all plotted particles
%   plotDataOut: structure containing plotted data.
%     Contains:
%	plotOutput.dataName ={'movieNo' 'particleNo' 'DDX' 'DDY' 't_Dex' 't_Aex' 'DD' 'DA' 'AD' 'AA' 'E' 'S'}
%	plotOutput.data - matrix containing columns specified in dataName. Each row is one frame
%   Ex, En - fret histogram data
%   Sx, Sn - stoiciometry histogram data 



n = numel(varargin);
i = 1;
nBin= 100;
eLim = [0 1];
sLim = [0 1];
filterArg = {};
selectParticleArg = {};
while i <= n
  if strcmp(varargin{i},'NBin')
    nBin = varargin{i+1};
    i=i+2;
  elseif strcmp(varargin{i},'ELim')
    eLim = varargin{i+1};
    i=i+2;
  elseif strcmp(varargin{i},'SLim')
    sLim = varargin{i+1};
    i=i+2;
  elseif iscell(varargin{i})
    filterArg = {filterArg{:}, varargin{i}};
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
sCol = find(strcmp(plotData.dataName,'S'));

ESmatrix = plotData.data(:,[eCol sCol]);

%calculate the limits for the plots
Min_E = eLim(1);
Max_E = eLim(2);
Min_S = sLim(1);
Max_S = sLim(2);

ESmatrix(find(ESmatrix(:,1)<Min_E | ESmatrix(:,1)>Max_E ...
		| ESmatrix(:,2)<Min_S | ESmatrix(:,2)>Max_S ), :) = []; 

%make sure we get exactly 100 bins by using bin centre
[eCentre eEdge]= bincentre(Min_E,Max_E,nBin);
[sCentre sEdge]= bincentre(Min_S,Max_S,nBin);

%plot it here
%-- Plot Bursts --
clf; fig1 = gcf; 
set(fig1,'color','w','name','Display Bursts')

subplot(3,3,[4 5 7 8]); fh1 = subplot(3,3,[4 5 7 8]); ax1 = gca;
Burst_Histogram = HIST2(ESmatrix(:,1), ESmatrix(:,2) ,eEdge,sEdge);
%hES = pcolor(eEdge, sEdge,Burst_Histogram); 
hES = pcolor(eCentre, sCentre,Burst_Histogram(2:end,2:end)); 
ax1 = gca;

load('colormap_gg.mat');
colormap(colormap_gg);

caxis([0, max(Burst_Histogram(:))]);

set(hES,'edgecolor','none');
%shading(ax1,'interp');
set(ax1,'TickDir', 'out');
xlabel('Epr','FontSize',14,'Parent',ax1);
ylabel('Sraw','FontSize',14,'Parent',ax1);

subplot(3,3,[1 2]); fh2 = subplot(3,3,[1 2]); ax2 = gca;
[TopN, TopX]= hist(ESmatrix(:,1),eCentre);
hE = bar(TopX,TopN);
Ex = TopX;
En = TopN;

ylabel('BurstX','FontSize',14,'Parent',ax2);
h01 = findobj(gca,'Type','patch');
set(h01,'FaceColor','b','EdgeColor','b');

subplot(3,3,[6 9]); fh3 =  subplot(3,3,[6 9]); ax3 = gca;
[SideN,SideX] = hist(ESmatrix(:,2) ,sCentre);
hS = barh(SideX,SideN);
Sx = SideX;
Sn = SideN;
h02 = findobj(gca,'Type','patch');
set(h02,'FaceColor','b','EdgeColor','b');
xlabel('BurstsX','FontSize',14,'Parent',ax3);

% Visual Touch-Ups
set(ax2,'XTick',[],'XTickLabel',[]);
set(ax3,'YTick',[],'YTickLabel',[]);

% set the axes labels and limits
axis(ax1,[Min_E Max_E Min_S Max_S]);
axis(ax2,[Min_E Max_E 0 max(TopN)*1.2]);
axis(ax3,[0 max(SideN)*1.2 Min_S Max_S]);

eTickLabels = [Min_E:((Max_E - Min_E)/5):Max_E];
sTickLabels = [Min_S:((Max_S - Min_S)/5):Max_S];
set(ax1,'XTick',eTickLabels, ...
	 'XTickLabel', eTickLabels, ...
     'YTick', sTickLabels, ...  
     'YTickLabel', sTickLabels);
set(ax1,'XTick',eTickLabels, ...
   	 'XTickLabel', eTickLabels, ...
  	 'YTick', sTickLabels, ...  
	 'YTickLabel', sTickLabels);
set(ax1,'XTick',eTickLabels, ...
   	 'XTickLabel', eTickLabels, ...
	 'YTick', sTickLabels, ...  
	 'YTickLabel', sTickLabels);

plotDataOut = plotData;


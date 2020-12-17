function [plotDataOut, absoluteIndex]= selectParticle(plotData,varargin)
%function [plotDataOut, absoluteIndex]= selectParticle(plotData,varargin)
% 
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 Alpha, released 101115
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%
n = numel(varargin);
i = 1;
if n >0
  useRelativeIndex = true;
  integrationTime = 1;
  while i <= n
    if strcmp(varargin{i},'RelativeIndex')
      useRelativeIndex = true;
      relativeIndex = varargin{i+1};
      i = i+2;
    elseif strcmp(varargin{i},'AbsoluteIndex')
      useRelativeIndex = false;
      absoluteIndex = varargin{i+1};
      i = i+2;
    else
      error('Unrecognised function argument');
    end
  end
  movieCol    = find(strcmp(plotData.dataName,'movieNo'));
  particleCol = find(strcmp(plotData.dataName,'particleNo'));

  % calculate the relative particle index
  if useRelativeIndex == true
    absoluteIndexList = unique(plotData.data(:,[movieCol,particleCol]),'rows');
    nParticle = size(absoluteIndexList,1);
    if relativeIndex <= nParticle
      absoluteIndex = absoluteIndexList(relativeIndex,:);
    else
      if nParticle == 1
        suffix = 'st';
      elseif nParticle == 2
        suffix = 'nd';
      elseif nParticle == 3
        suffix = 'rd';
      else
        suffix = 'th';
      end

      error(['Attempted to access ',num2str(relativeIndex),suffix, ' particle, but there are only ',num2str(nParticle), ' particles.']);
    end
  end
  movieIndex = absoluteIndex(:,1);
  particleIndex = absoluteIndex(:,2);

  % filter selected particles
  matchRow = find( ismember(plotData.data(:,movieCol),movieIndex) & ismember(plotData.data(:,particleCol),particleIndex));
  plotDataOut = plotData;
  plotDataOut.data = plotDataOut.data(matchRow,:);
else% return all rows
  movieCol    = find(strcmp(plotData.dataName,'movieNo'));
  particleCol = find(strcmp(plotData.dataName,'particleNo'));
  absoluteIndex =  unique(plotData.data(:,[movieCol,particleCol]),'rows');
  plotDataOut = plotData;
end

if isempty(plotDataOut.data)
  error('No data meets these criteria')
end

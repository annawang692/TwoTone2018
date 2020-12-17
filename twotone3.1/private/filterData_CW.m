function plotOutput  = filterData_CW(file, varargin)
% function plotOutput  = filterData_CW(file, varargin)
% INPUT:
%   file: either a file filter, eg '*.awesomedata.mat' or a cell of filenames {'awesomedata1.mat', 'awesomedate2.mat, ...}
%   Optional arguments:
%     filter (any number of these may be applied): {'filterType','filterParam','>' OR '<', Value}
%     'RelativeParticleIndex', vector, eg, [1, 4, 6] to plot only a certain # of the molecules
%	filter: {FilterType, FilterParam, Range, Value}
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
% OUTPUT 
%   movieNo particleNo DX DY t D A E
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
i =1;
k = 1;
selectParticle = false;
filter = {};
while i <= n
  if iscell(varargin{i})
    filter{k} = varargin{i};
    k = k + 1;
    i = i+1;
  end
end

twotoneData = loadData(file); %returns a 1xn structure containing each file
nFilter = numel(filter);
% for each particle in each movie, return a list of all frames satisfying the filter
frameList = applyFilter(twotoneData,filter);
% given a list of all frames, produce a vector of output data for each particle
fretDataName = twotoneData.settings.fretDataName;
plotOutput.dataName = {'movieNo',  'particleNo' 'DX' 'DY', fretDataName{:} };
plotOutput.data = convertDataToMat(twotoneData,frameList);

%---------------------------------------------------------------
function matchrowsparticle = filterParticle(data,filter, fretDataName, aDetParamName, dcol,acol)
% filter: {FilterType, FilterParam, Range, Value}
% FilterType:  'fretData', 
%   FilterParam:   't','D','A','E'
% FilterType: 'all_aDetParam'
%   FilterParam:  'nearestNeighbor'
% FilterType:    'D_aDetParam'
% FilterType:    'A_aDetParam'
%   FilterParam:     'amplitude', 'sx', 'sy', 'eccentricity', 
%   
% FilterType:    'D_aDetPos'
% FilterType:    'A_aDetPos'
%   FilterParam:     'X' 'Y'
% Range: '<' or '>'
% Value: Scalar number
%
% match filter data

allrowsparticle = find(data.fretData(:,1));
if strcmp(filter{1}, 'fretData')
  fretFilterCol = find(strcmp(fretDataName, filter{2}));
  if strcmp(filter{3},'<') 
      matchrowsparticle = find(data.fretData(:,fretFilterCol) < filter{4});
  elseif strcmp(filter{3},'>') 
      matchrowsparticle = find(data.fretData(:,fretFilterCol) > filter{4});
  else
      matchrowsparticle=[];
  end
elseif strcmp(filter{1}, 'all_aDetParam')
  if strcmp(filter{2},'nearestNeighbor')
    if strcmp(filter{3},'<') && (data.aDetData.nearestNeighborDist < filter{4})
      matchrowsparticle = allrowsparticle;
    elseif strcmp(filter{3},'>') && (data.aDetData.nearestNeighborDist > filter{4})
      matchrowsparticle = allrowsparticle;
    else
      matchrowsparticle=[];
    end
  end
elseif strcmp(filter{1}, 'D_aDetParam')||strcmp(filter{1}, 'A_aDetParam')
  if strcmp(filter{1}, 'D_aDetParam')
    aDetRow = dcol;
  elseif strcmp(filter{1}, 'A_aDetParam')
    aDetRow = acol;
  end

  aDetFilterCol = find(strcmp(aDetParamName, filter{2}));
  if strcmp(filter{3},'<') && (data.aDetData.aDetParam(aDetRow,aDetFilterCol) < filter{4})
      matchrowsparticle = allrowsparticle;
  elseif strcmp(filter{3},'>')  && (data.aDetData.aDetParam(aDetRow,aDetFilterCol) < filter{4})
      matchrowsparticle = allrowsparticle;
  else
      matchrowsparticle=[];
  end
elseif strcmp(filter{1}, 'D_aDetPos')||strcmp(filter{1}, 'A_aDetPos')
  if strcmp(filter{1}, 'D_aDetPos')
    aDetRow = dcol;
  elseif strcmp(filter{1}, 'A_aDetPos')
    aDetRow = acol;
  end
  
  if strcmp(filter{2}, 'X')
    posCol = 1;
  elseif strcmp(filter{2}, 'Y')
    posCol = 2;
  end
  
  if strcmp(filter{3},'<') && (data.aDetData.aDetParam(aDetRow,posCol) < filter{4})
      matchrowsparticle = allrowsparticle;
  elseif strcmp(filter{3},'>')  && (data.aDetData.aDetParam(aDetRow,posCol) < filter{4})
      matchrowsparticle = allrowsparticle;
  else
      matchrowsparticle=[];
  end
else
  error('Unrecognised filter option');
end

%----------------------------------------------------------------------------
function twotoneData = loadData(file)

fileCell = parseFileFilter(file);

for i = 1:numel(fileCell)
  temp = load(fileCell{i},'twotoneData');
  twotoneData(i) = temp.twotoneData;
end

%-----------------------------------------------------------------------------
function frameList = applyFilter(twotoneData,filter);
% function frameList = applyFilter(twotoneData,filter);
%  for each particle in each movie, return a list of all frames satisfying the filter

fretDataName = twotoneData.settings.fretDataName;
aDetParamName = twotoneData.results.aDetParamName;
dcol = strcmp(twotoneData.settings.imageSettings.aDetChannelName,'D');
acol = strcmp(twotoneData.settings.imageSettings.aDetChannelName,'A');

% TEMP findallrows
findallrows = @(x) find(x.fretData(:,1));
findallrowseachmovie = @(x)  arrayfun(@(y) (findallrows(y)), x.results.data,'UniformOutput',false);
findallrowsallmovie = @(x) arrayfun(@(y) (findallrowseachmovie(y)), x,'UniformOutput',false);

findfilteredrowseachmovie = @(x,filter, fretDataName, aDetParamName, dcol,acol) ...
    arrayfun(@(y)  (filterParticle(y,filter, fretDataName, aDetParamName, dcol,acol)), x.results.data,'UniformOutput',false);
findfilteredrowsallmovie = @(x,filter, fretDataName, aDetParamName, dcol,acol) ...
    arrayfun(@(y) (findfilteredrowseachmovie(y,filter, fretDataName, aDetParamName, dcol,acol)), x,'UniformOutput',false);

%initial list returns everything
frameList= findallrowsallmovie(twotoneData);

intersecttwontoneparticle = @(list1,list2) cellfun(@(list1,list2) intersect(list1,list2), list1,list2, 'UniformOutput',false);
intersecttwotonemovie = @(list1,list2) cellfun(@(list1,list2) intersecttwontoneparticle(list1,list2), list1,list2, 'UniformOutput',false);
nFilter = numel(filter);
for i=1:nFilter
  filterMatchRows = findfilteredrowsallmovie(twotoneData,filter{i}, fretDataName, aDetParamName, dcol,acol);
  frameList = intersecttwotonemovie(frameList,filterMatchRows);
end

%-----------------------------------------------------------------------------
function plotOutput = convertDataToMat(twotoneData,frameList)
% given list of all frames, produce an array of output data for each particle
% output is
%   movieNo particleNo DX DY t D A E 
frameListStruct = cell2struct(frameList(:)','data',1);

nMovie = numel(frameList);
movieNo = 1:nMovie;
plotOutCell = arrayfun(@(x,frameListStruct,movieNo) makePlotArrayMovie(x.results.data,frameListStruct,movieNo), ...
		twotoneData(:),frameListStruct(:),movieNo(:),'UniformOutput',false);
plotOutput = cell2mat(plotOutCell);



%-----------------------------------------------------------------------------
function plotOut = makePlotArrayParticle(dataParticle,frameList,particleNo,movieNo, relParticleStart)
% function to make an array for each particle in a movie

nRow = size(dataParticle.fretData,2)+4;
if ~isempty(frameList.data)
  nCol = numel(frameList.data);
  plotOut = zeros(nCol,nRow);
  plotOut(:,1) = movieNo;
  plotOut(:,2) = particleNo;
  DposX =dataParticle.aDetData.aDetPos(1,1);
  DposY =dataParticle.aDetData.aDetPos(1,2);
  plotOut(:,3) = DposX;
  plotOut(:,4) = DposY;
  plotOut(:,5:end) = dataParticle.fretData(frameList.data,:);
else
  plotOut = zeros(0,nRow);
end

%-----------------------------------------------------------------------------
function plotOut = makePlotArrayMovie(dataMovie,frameList,movieNo)
%function to make an array for all particles in a movie
% @(x,filter) arrayfun(@(y) (findallrowseachmovie(y,filter)), x,'UniformOutput',false);

nParticle = numel(dataMovie);
particleNo = (1:nParticle)';
frameListStruct = cell2struct(frameList.data(:)','data',1);
plotOutCell = arrayfun(@(y,frameListStruct,particleNo) (makePlotArrayParticle(y,frameListStruct,particleNo,movieNo)), dataMovie,frameListStruct,particleNo,'UniformOutput',false);
plotOut = cell2mat(plotOutCell);

























function plotOutput  = filterData_ALEX(file, varargin)
% function plotOutput  = filterData_ALEX(file, varargin)
% INPUT:
%   file: either a file filter, eg '*.awesomedata.mat' or a cell of filenames {'awesomedata1.mat', 'awesomedate2.mat, ...}
%   Optional arguments:
%     filter (any number of these may be applied): {'filterType','filterParam','>' OR '<', Value}
%	filter: {FilterType, FilterParam, Range, Value}
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
% OUTPUT 
%   movieNo particleNo DDX DDY t_Dex t_Aex DD DA AD AA E S
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
 % elseif strcmp(varargin{i},'RelativeParticleIndex')
 %   selectParticle = true;
 %   particleVector = varargin{i+1};
 %   i = i +2;
  else
    varargin{i}
    error('Unrecognised filter argument');
  end
end

twotoneData = loadData(file); %returns a 1xn structure containing each file
nFilter = numel(filter);
% for each particle in each movie, return a list of all frames satisfying the filter
frameList = applyFilter(twotoneData,filter);
% given a list of all frames, produce a vector of output data for each particle
fretDataName = twotoneData.settings.fretDataName;
plotOutput.dataName = {'movieNo',  'particleNo' 'DDX' 'DDY', fretDataName{:} };
plotOutput.data = convertDataToMat(twotoneData,frameList);

%---------------------------------------------------------------
function matchrowsparticle = filterParticle(data,filter, fretDataName, aDetParamName, ddcol,aacol,dacol)
% filter: {FilterType, FilterParam, Range, Value}
% FilterType:  'fretData', 
%   FilterParam:   't_Dex', 't_Aex', 'DD', 'DA','AA', 'AA', 'E', 'S'
% FilterType: 'all_aDetParam'
%   FilterParam:  'nearestNeighbor'
% FilterType:    'DD_aDetParam'
% FilterType:    'DA_aDetParam'
% FilterType:    'AA_aDetParam'
%   FilterParam:     'amplitude', 'sx', 'sy', 'eccentricity', 
%   
% FilterType:    'DD_aDetPos'
% FilterType:    'DA_aDetPos'
% FilterType:    'AA_aDetPos'
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
elseif strcmp(filter{1}, 'DD_aDetParam')||strcmp(filter{1}, 'AA_aDetParam')||strcmp(filter{1}, 'DA_aDetParam')
  if strcmp(filter{1}, 'DD_aDetParam')
    aDetRow = ddcol;
  elseif strcmp(filter{1}, 'AA_aDetParam')
    aDetRow = aacol;
  elseif strcmp(filter{1}, 'DA_aDetParam')
    aDetRow = dacol;
  end

  aDetFilterCol = find(strcmp(aDetParamName, filter{2}));
  if strcmp(filter{3},'<') && (data.aDetData.aDetParam(aDetRow,aDetFilterCol) < filter{4})
      matchrowsparticle = allrowsparticle;
  elseif strcmp(filter{3},'>')  && (data.aDetData.aDetParam(aDetRow,aDetFilterCol) < filter{4})
      matchrowsparticle = allrowsparticle;
  else
      matchrowsparticle=[];
  end
elseif strcmp(filter{1}, 'DD_aDetPos')||strcmp(filter{1}, 'AA_aDetPos')||strcmp(filter{1}, 'DA_aDetPos')
  if strcmp(filter{1}, 'DD_aDetPos')
    aDetRow = ddcol;
  elseif strcmp(filter{1}, 'AA_aDetPos')
    aDetRow = aacol;
  elseif strcmp(filter{1}, 'DA_aDetPos')
    aDetRow = dacol;
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
ddcol = strcmp(twotoneData.settings.imageSettings.aDetChannelName,'DexDem');
aacol = strcmp(twotoneData.settings.imageSettings.aDetChannelName,'AexAem');
dacol = strcmp(twotoneData.settings.imageSettings.aDetChannelName,'DexAem');

% TEMP findallrows
findallrows = @(x) find(x.fretData(:,1));
findallrowseachmovie = @(x)  arrayfun(@(y) (findallrows(y)), x.results.data,'UniformOutput',false);
findallrowsallmovie = @(x) arrayfun(@(y) (findallrowseachmovie(y)), x,'UniformOutput',false);

findfilteredrowseachmovie = @(x,filter, fretDataName, aDetParamName, ddcol,aacol,dacol) ...
    arrayfun(@(y)  (filterParticle(y,filter, fretDataName, aDetParamName, ddcol,aacol,dacol)), x.results.data,'UniformOutput',false);
findfilteredrowsallmovie = @(x,filter, fretDataName, aDetParamName, ddcol,aacol,dacol) ...
    arrayfun(@(y) (findfilteredrowseachmovie(y,filter, fretDataName, aDetParamName, ddcol,aacol,dacol)), x,'UniformOutput',false);

%initial list returns everything
frameList= findallrowsallmovie(twotoneData);

intersecttwontoneparticle = @(list1,list2) cellfun(@(list1,list2) intersect(list1,list2), list1,list2, 'UniformOutput',false);
intersecttwotonemovie = @(list1,list2) cellfun(@(list1,list2) intersecttwontoneparticle(list1,list2), list1,list2, 'UniformOutput',false);
nFilter = numel(filter);
for i=1:nFilter
  filterMatchRows = findfilteredrowsallmovie(twotoneData,filter{i}, fretDataName, aDetParamName, ddcol,aacol,dacol);
  frameList = intersecttwotonemovie(frameList,filterMatchRows);
end

%-----------------------------------------------------------------------------
function plotOutput = convertDataToMat(twotoneData,frameList)
% given list of all frames, produce an array of output data for each particle
% output is
%   movieNo particleNo DDX DDY t_Dex t_Aex DD DA AD AA E S
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
  DDposX =dataParticle.aDetData.aDetPos(1,1);
  DDposY =dataParticle.aDetData.aDetPos(1,2);
  plotOut(:,3) = DDposX;
  plotOut(:,4) = DDposY;
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



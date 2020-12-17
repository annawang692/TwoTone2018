function [filteredClusters, allValidClusters,distanceDistribution] = associateCWchannels(D,A, xmin,Dparams,Aparams,paramNames,filterChoice, TFORM)
%function [filteredClusters, allValidClusters,distanceDistribution] = associateCWchannels(D,A, xmin,Dparams,Aparams,paramNames,filterChoice, TFORM)
% data is clustered based on distance and then filtered based on stoichtiometry. ambiguous molecules are excluded (>1 particle in any channel)
%
% filteredClusters:
%         NNcur: 0.0100
%           nd: 0
%           na: 0
%            d: []
%            a: []
%      dparams: []
%      aparams: []
%    paramNames: {}
%
% distanceDistribution is the cluster distance distribution
%
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 Alpha, released 101115
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%

if ~exist('TFORM')  
  TFORM =[];
end

tformIsGreenToRed = 1;
[D,A] = tformChannels(D, A,TFORM, tformIsGreenToRed); %transform the red positions into the green coordinate system

[validClusters, distanceDistribution  ] = clusterChannels(D,A, xmin,Dparams,Aparams,paramNames);%cluster the data based on distance, reject molecules with inappropriate stoichiomety

[validClusters] = invtformClusters(validClusters,TFORM,tformIsGreenToRed);%transform the red positions back into the red coordinate system
allValidClusters = validClusters;

[filteredClusters] = filterChannels(validClusters, filterChoice); %filter based on stoichiomety


%-------------------------------------------------------------------------------------------------------
function [validClusters, distanceDistribution  ] = clusterChannels(D,A, xmin,Dparams,Aparams,paramNames)
%function [validClusters, distanceDistribution  ] = clusterChannels(D,A, xmin,Dparams,Aparams,paramNames)
% group the channels based on distance 
% stoichiometries >1 in any channel are rejected 
% validClusters contains: NN, 
%		d = [x,y],a = [x,y]  - points within each cluster
%		 n_d, n_a,  - number of each point
% distance distribution is the linkage output


%1.form set of all points X = {D,A}
nD = size(D,1);
nA = size(A,1);
nPoints = nD+nA;
X = [D; A];

%sort out the params data
if ~exist('Dparams') || isempty(Dparams)
  Dparams = zeros(nD,0);%empty matrix which wont cause crashes
end  
if ~exist('Aparams') || isempty(Aparams)
  Aparams = zeros(nA,0);%empty matrix which wont cause crashes
end  
Xparams = [Dparams;Aparams];

%2. Form clusters C
% cluster the data based on distance
Y = pdist(X,'euclid'); 
Z = linkage(Y,'single'); 
distanceDistribution = Z(:,3);
T = cluster(Z,'cutoff',xmin,'criterion','distance');%T is the number of the cluster to which the particle belongs

NNout = particleNN(Z,xmin);% NNout is the nearest neighbour distance of each cluster
channelVector = [1.*ones(nD,1); 2.*ones(nA,1)]; % 1 means D, 2 means A
Xclust = [X,T,NNout, channelVector]; %append the cluster number & nearest neighbour distance and which channel its from
% Xclust : [x, y , clNum, NNout, channel]

%3. For each Ci, calculate nearest neighbour distance NN, and number of members n_d, n_a, 
% Ci contains: NN, 
%		d = [x,y],a = [x,y] - points within each cluster
%		 n_d, n_a, - number of each point
nClusters = max(T);
for i = 1:nClusters
  clRows = find(Xclust(:,3) == i);
  C(i).NNcur = Xclust(clRows(1),4);
  C(i).nd = 0;
  C(i).na = 0;
  C(i).d=zeros(0,2);
  C(i).a=zeros(0,2);
  C(i).dparams=zeros(0,4);
  C(i).aparams=zeros(0,4);
  C(i).paramNames=paramNames;
  for j =1 :numel(clRows)
    curRow = clRows(j);
    if Xclust(curRow,5) == 1
      C(i).d = [C(i).d; Xclust(curRow,1:2)];
      C(i).dparams=[C(i).dparams;Xparams(curRow,:)];
      C(i).nd = C(i).nd + 1;
    elseif Xclust(curRow,5) == 2
      C(i).a = [C(i).a; Xclust(curRow,1:2)];
      C(i).aparams=[C(i).aparams;Xparams(curRow,:)];
      C(i).na = C(i).na + 1;
    end
  end
end 

%4. Form sets validClusters   ={Ci(n_d<=1 & n_a<=1 )}
%	      invalidClusters ={Ci(n_d>1  & n_a>1 )}
k = 1;
for i = 1:nClusters
  if ((C(i).nd <= 1) &&  (C(i).na <= 1) )
    validClusters(k) = C(i);
    k = k+1;
  end
end

nParam = numel(validClusters(1).paramNames);
%ad the positions in the empty cluster subfields
nValidCluster = numel(validClusters);
for i = 1:nValidCluster
  if isempty(validClusters(i).d)
    validClusters(i).dparams=NaN*ones(1,nParam);
    validClusters(i).d= validClusters(i).a;
  elseif isempty(validClusters(i).a)
    validClusters(i).aparams=NaN*ones(1,nParam);
    validClusters(i).a= validClusters(i).d; 
  end
end

%----------------------------------------------------
function NNout = particleNN(Z,xMin)
%function NNout = particleNN(Z,xMin)
%
% return the nearest neighbour distances for the clustered particles (type of NN returned depends on 
% linkage method used to form Z
% Z is linkage output, and xMin is the cluster cutoff

maxDistance = max(Z(:,3));
nPoints = size(Z,1)+1;
if xMin > maxDistance%ie only one cluster found
  NNout = ones(nPoints,1).*Inf;
  warning('xMin exceeded largest nearest neighbour distance');
else
  crit = Z(:,3);
  conn = find(crit>xMin);%above theshold-distance clusters

  tempDistances=[];

  for i = 1:numel(conn) 
    currentDistance = Z(conn(i),3);
    for j = 1:2 %left and right columns
      currentCluster = Z(conn(i),j);
      allChildren = findChildren(Z, currentCluster);%all of the children (cluster members) have this distance too
      nChildren = numel(allChildren); 
      tempDistances = [tempDistances; allChildren , ones(nChildren,1).*currentDistance];%first above-thresh NN distance 
    end
  end

  %now that we have a list of distances, we need to remove the repeats
  for i = 1:nPoints
    allDistanceVal = tempDistances(find(tempDistances(:,1)==i) , 2); %get all the nnVal for point i 
    nnVal = min(allDistanceVal);
    NNout(i,1) = nnVal;
  end
end

%----------------------
function a= findChildren(Z,b)
%b is parent cluster number
%return all child clusters

nPoints = size(Z,1)+1;
if b > (2*nPoints -1)
  error('b is larger than the total number of clusters');
else
  if b <=nPoints
    a = b;
  else
    a = [];
    immediateChildren = Z(b-nPoints,1:2);
    for j=1:2%left and right columns
      cTemp = immediateChildren(j);
      if cTemp <=nPoints
	a = [a; cTemp];
      else
	cOut = findChildren(Z,cTemp);
	a = [a; cOut];
      end
    end
  end
end

%---------------------------------------------------------------------
function [Dtform, Atform] = tformChannels(D,A,TFORM,tformIsGreenToRed);

if ~isempty(TFORM)% we have a valid TFORM
  if tformIsGreenToRed == 1%tform green channel
    Dtform = tformfwd(TFORM,D);
    Atform = A;
  else	%tform red channel
    Dtform = D;
    Atform = tformfwd(TFORM,A);
  end
else%no tform supplied so do nothing
  Dtform = D;
  Atform = A;
end

%---------------------------------------------------------------------
function [cOut] = invtformClusters(cIn,TFORM,tformIsGreenToRed);

cOut = cIn;

nCluster = numel(cOut);
if ~isempty(TFORM)
  for i = 1:nCluster
    if tformIsGreenToRed == 1%tform green channel
      cOut(i).d = tforminv(TFORM,cOut(i).d);
      cOut(i).a = cOut(i).a;
    else %tform red channel
      cOut(i).d = cOut(i).d;
      cOut(i).a = tforminv(TFORM,cOut(i).a);
    end
  end
else
  %no tform supplied so do nothing
end

%---------------------------------------------------------------------
function [filteredClusters] = filterChannels(validClusters, filterChoice); %filter based on stoichiomety

nCluster = numel(validClusters);
k=1;
for i = 1:nCluster
  
  %decide whether or not to include the point based on the filter supplied
  if strcmp(filterChoice,'showAll')
    adThisCluster = true;
  elseif strcmp(filterChoice,'D') && validClusters(i).nd==1
    adThisCluster = true;
  elseif strcmp(filterChoice,'A') && validClusters(i).na==1
    adThisCluster = true;
  elseif strcmp(filterChoice,'D&&A') && validClusters(i).nd==1 && validClusters(i).na==1
    adThisCluster = true;
  elseif strcmp(filterChoice,'D||A') && (validClusters(i).nd==1 || validClusters(i).na==1)
    adThisCluster = true;
  else
    adThisCluster = false;
  end
  
  if adThisCluster == true
    filteredClusters(k) = validClusters(i);
    k = k + 1;
  end
end

if ~exist('filteredClusters')
  filteredClusters =[];
end

function [filteredClusters, allValidClusters,distanceDistribution] = associateALEXchannels(DD,AA,DA, xmin,DDparams,AAparams,DAparams,paramNames,filterChoice, TFORM)
%function [filteredClusters, allValidClusters,distanceDistribution] = associateALEXchannels(DD,AA,DA, xmin,DDparams,AAparams,DAparams,paramNames,filterChoice, TFORM)
% data is clustered based on distance and then filtered based on stoichtiometry. ambiguous molecules are excluded (>1 particle in any channel)
%
% filteredClusters:
%         NNcur: 0.0100
%           ndd: 0
%           naa: 0
%          nda: 1
%            dd: []
%            aa: []
%            da: [0.9716 0.6150]
%      ddparams: []
%      aaparams: []
%      daparams: [1x0 double]
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
[DD,AA,DA] = tformChannels(DD, AA,DA,TFORM, tformIsGreenToRed); %transform the red positions into the green coordinate system

[validClusters, distanceDistribution  ] = clusterChannels(DD,AA,DA, xmin,DDparams,AAparams,DAparams,paramNames);%cluster the data based on distance, reject molecules with inappropriate stoichiomety

[validClusters] = invtformClusters(validClusters,TFORM,tformIsGreenToRed);%transform the red positions back into the red coordinate system
allValidClusters = validClusters;

[filteredClusters] = filterChannels(validClusters, filterChoice); %filter based on stoichiomety


%-------------------------------------------------------------------------------------------------------
function [validClusters, distanceDistribution  ] = clusterChannels(DD,AA,DA, xmin,DDparams,AAparams,DAparams,paramNames)
%function [validClusters, invalidClusters,distanceDistribution  ] = clusterChannels(DD,AA,DA, xmin,DDparams,AAparams,DAparams)
% group the channels based on distance 
% stoichiometries >1 in any channel are rejected 
% validClusters contains: NN, 
%		dd = [x,y],aa = [x,y] ,da = [x,y] - points within each cluster
%		 n_dd, n_aa, n_da - number of each point
% distance distribution is the linkage output


%1.form set of all points X = {DD,AA,DA}
nDD = size(DD,1);
nAA = size(AA,1);
nDA = size(DA,1);
nPoints = nDD+nAA+nDA;
X = [DD; AA ; DA];

%sort out the params data
if ~exist('DDparams') || isempty(DDparams)
  DDparams = zeros(nDD,0);%empty matrix which wont cause crashes
end  
if ~exist('AAparams') || isempty(AAparams)
  AAparams = zeros(nAA,0);%empty matrix which wont cause crashes
end  
if ~exist('DAparams') || isempty(DAparams)
  DAparams = zeros(nDA,0);%empty matrix which wont cause crashes
end  
Xparams = [DDparams;AAparams;DAparams];

%2. Form clusters C
% cluster the data based on distance
Y = pdist(X,'euclid'); 
Z = linkage(Y,'single'); 
distanceDistribution = Z(:,3);
T = cluster(Z,'cutoff',xmin,'criterion','distance');%T is the number of the cluster to which the particle belongs

NNout = particleNN(Z,xmin);% NNout is the nearest neighbour distance of each cluster
channelVector = [1.*ones(nDD,1); 2.*ones(nAA,1); 3*ones(nDA,1)]; % 1 means DD, 2 means AA, 3 means DA
Xclust = [X,T,NNout, channelVector]; %append the cluster number & nearest neighbour distance and which channel its from
% Xclust : [x, y , clNum, NNout, channel]

%3. For each Ci, calculate nearest neighbour distance NN, and number of members n_dd, n_aa, n_da
% Ci contains: NN, 
%		dd = [x,y],aa = [x,y] ,da = [x,y] - points within each cluster
%		 n_dd, n_aa, n_da - number of each point
nClusters = max(T);
for i = 1:nClusters
  clRows = find(Xclust(:,3) == i);
  C(i).NNcur = Xclust(clRows(1),4);
  C(i).ndd = 0;
  C(i).naa = 0;
  C(i).nda = 0;
  C(i).dd=zeros(0,2);
  C(i).aa=zeros(0,2);
  C(i).da=zeros(0,2);
  C(i).ddparams=zeros(0,4);
  C(i).aaparams=zeros(0,4);
  C(i).daparams=zeros(0,4);
  C(i).paramNames=paramNames;
  for j =1 :numel(clRows)
    curRow = clRows(j);
    if Xclust(curRow,5) == 1
      C(i).dd = [C(i).dd; Xclust(curRow,1:2)];
      C(i).ddparams=[C(i).ddparams;Xparams(curRow,:)];
      C(i).ndd = C(i).ndd + 1;
    elseif Xclust(curRow,5) == 2
      C(i).aa = [C(i).aa; Xclust(curRow,1:2)];
      C(i).aaparams=[C(i).aaparams;Xparams(curRow,:)];
      C(i).naa = C(i).naa + 1;
    elseif Xclust(curRow,5) == 3
      C(i).da = [C(i).da; Xclust(curRow,1:2)];
      C(i).daparams=[C(i).daparams;Xparams(curRow,:)];
      C(i).nda = C(i).nda + 1;
    end
  end
end 

%4. Form sets validClusters   ={Ci(n_dd<=1 & n_aa<=1 & n_da<=1)}
%	      invalidClusters ={Ci(n_dd>1  & n_aa>1  & n_da>1)}
k = 1;
for i = 1:nClusters
  if (C(i).ndd <= 1) &&  (C(i).naa <= 1) && (C(i).nda <= 1)
    validClusters(k) = C(i);
    k = k+1;
  end
end

nParam = numel(validClusters(1).paramNames);
%add the positions in the empty cluster subfields
nValidCluster = numel(validClusters);
for i = 1:nValidCluster
  if isempty(validClusters(i).dd)
    validClusters(i).ddparams=NaN*ones(1,nParam);
    if ~isempty(validClusters(i).aa) && ~isempty(validClusters(i).da)
      ddFinal = mean( [validClusters(i).aa; validClusters(i).da],1);
    elseif isempty(validClusters(i).aa) && ~isempty(validClusters(i).da)
      ddFinal = validClusters(i).da;
    elseif ~isempty(validClusters(i).aa) && isempty(validClusters(i).da)
      ddFinal = validClusters(i).aa;
    end
  else 
    ddFinal = validClusters(i).dd;
  end

  if isempty(validClusters(i).aa)
    validClusters(i).aaparams=NaN*ones(1,nParam);
    if ~isempty(validClusters(i).dd) && ~isempty(validClusters(i).da)
      aaFinal = validClusters(i).da; 
    elseif isempty(validClusters(i).dd) && ~isempty(validClusters(i).da)
      aaFinal = validClusters(i).da;
    elseif ~isempty(validClusters(i).dd) && isempty(validClusters(i).da)
      aaFinal = validClusters(i).dd;
    end
  else 
    aaFinal = validClusters(i).aa;
  end

  if isempty(validClusters(i).da)
    validClusters(i).daparams=NaN*ones(1,nParam);
    if ~isempty(validClusters(i).dd) && ~isempty(validClusters(i).aa)
      daFinal = validClusters(i).aa; 
    elseif isempty(validClusters(i).dd) && ~isempty(validClusters(i).aa)
      daFinal = validClusters(i).aa;
    elseif ~isempty(validClusters(i).dd) && isempty(validClusters(i).aa)
      daFinal = validClusters(i).dd;
    end
  else 
    daFinal = validClusters(i).da;
  end
  
  validClusters(i).dd = ddFinal;
  validClusters(i).aa = aaFinal;
  validClusters(i).da = daFinal;
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
function [DDtform, AAtform,DAtform] = tformChannels(DD,AA,DA,TFORM,tformIsGreenToRed);

if ~isempty(TFORM)% we have a valid TFORM
  if tformIsGreenToRed == 1%tform green channel
    DDtform = tformfwd(TFORM,DD);
    AAtform = AA;
    DAtform = DA;
  else	%tform red channel
    DDtform = DD;
    AAtform = tformfwd(TFORM,AA);
    DAtform = tformfwd(TFORM,DA);
  end
else%no tform supplied so do nothing
  DDtform = DD;
  AAtform = AA;
  DAtform = DA;
end

%---------------------------------------------------------------------
function [cOut] = invtformClusters(cIn,TFORM,tformIsGreenToRed);

cOut = cIn;

nCluster = numel(cOut);
if ~isempty(TFORM)
  for i = 1:nCluster
    if tformIsGreenToRed == 1%tform green channel
      cOut(i).dd = tforminv(TFORM,cOut(i).dd);
      cOut(i).aa = cOut(i).aa;
      cOut(i).da = cOut(i).da;
    else %tform red channel
      cOut(i).dd = cOut(i).dd;
      cOut(i).aa = tforminv(TFORM,cOut(i).aa);
      cOut(i).da = tforminv(TFORM,cOut(i).da);
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
    addThisCluster = true;
  elseif strcmp(filterChoice,'DexDem') && validClusters(i).ndd==1
    addThisCluster = true;
  elseif strcmp(filterChoice,'AexAem') && validClusters(i).naa==1
    addThisCluster = true;
  elseif strcmp(filterChoice,'DexAem') && validClusters(i).nda==1
    addThisCluster = true;
  elseif strcmp(filterChoice,'DexDem&&AexAem') && validClusters(i).ndd==1 && validClusters(i).naa==1
    addThisCluster = true;
  elseif strcmp(filterChoice,'DexDem&&DexAem') && validClusters(i).ndd==1 && validClusters(i).nda==1
    addThisCluster = true;
  elseif strcmp(filterChoice,'AexAem&&DexAem') && validClusters(i).naa==1 && validClusters(i).nda==1
    addThisCluster = true;
  elseif strcmp(filterChoice,'DexDem&&AexAem&&DexAem') && validClusters(i).ndd==1 && validClusters(i).naa==1 && validClusters(i).nda==1
    addThisCluster = true;
  elseif strcmp(filterChoice,'DexDem||AexAem||DexAem') && (validClusters(i).ndd==1 || validClusters(i).naa==1 || validClusters(i).nda==1)
    addThisCluster = true;
  elseif strcmp(filterChoice,'DexDem||AexAem') && (validClusters(i).ndd==1 || validClusters(i).naa==1)
    addThisCluster = true;
  else
    addThisCluster = false;
  end
  
  if addThisCluster == true
    filteredClusters(k) = validClusters(i);
    k = k + 1;
  end
end

if ~exist('filteredClusters')
  filteredClusters =[];
end

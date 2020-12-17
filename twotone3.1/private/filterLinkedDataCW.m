function [filteredCluster, Dlinkedpos, Alinkedpos] = filterLinkedDataCW(clusterIn, params);
%function [filteredCluster, Dlinkedpos, Alinkedpos] = filterLinkedDataCW(clusterIn, params);
% filter the data based on nearest neighbour distance etc
%params.applyNN		
%params.applyEccentricity	
%params.applySigma		
%params.NNlim			
%params.eccLim		
%params.sigmaLim	
% 
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 Alpha, released 101115
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%

nCluster = numel(clusterIn);
if nCluster > 0  
  eccCol = find(strcmp(clusterIn(1).paramNames,'eccentricity'));
  sColX  = find(strcmp(clusterIn(1).paramNames,'sx'));
  sColY= find(strcmp(clusterIn(1).paramNames,'sy'));
  k =1;
  for i = 1:nCluster
    currCluster = clusterIn(i);
    nnFail = (params.applyNN == true && (currCluster.NNcur < params.NNlim) );
    
    eccFail = (params.applyEccentricity== true ...
	&&(     (currCluster.dparams(eccCol) > params.eccLim(1))  ...
	    ||  (currCluster.aparams(eccCol) > params.eccLim(1)) ) );
    sigmaXFail = (params.applySigma== true ...
	&&(    (currCluster.dparams(sColX) < params.sigmaLim(1)) ||(currCluster.dparams(sColX) > params.sigmaLim(2)) ...
	    || (currCluster.aparams(sColX) < params.sigmaLim(1)) ||(currCluster.aparams(sColX) > params.sigmaLim(2)) ) );
    sigmaYFail = (params.applySigma== true ...
	&&(    (currCluster.dparams(sColY) < params.sigmaLim(1)) ||(currCluster.dparams(sColY) > params.sigmaLim(2)) ...
	    || (currCluster.aparams(sColY) < params.sigmaLim(1)) ||(currCluster.aparams(sColY) > params.sigmaLim(2)) ) );
    if ~( nnFail || eccFail || sigmaXFail || sigmaYFail)
      filteredCluster(k) = currCluster;
      k=k+1;
    end
    
  end
else
  filteredCluster = clusterIn;
end

if ~exist('filteredCluster')
  filteredCluster =[];
end
Dlinkedpos=zeros(0,2);
Alinkedpos=zeros(0,2);

nClusterOut = numel(filteredCluster);
for i = 1:nClusterOut
  Dlinkedpos=[Dlinkedpos;filteredCluster(i).d];
  Alinkedpos=[Alinkedpos;filteredCluster(i).a];
end



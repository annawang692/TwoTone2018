function [filteredCluster, DDlinkedpos, AAlinkedpos,DAlinkedpos] = filterLinkedDataALEX(clusterIn, params);
%function [filteredCluster, DDlinkedpos, AAlinkedpos,DAlinkedpos] = filterLinkedDataALEX(clusterIn, params);
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
	&&(     (currCluster.ddparams(eccCol) > params.eccLim(1))  ...
	    ||  (currCluster.aaparams(eccCol) > params.eccLim(1))  ...
	    ||  (currCluster.daparams(eccCol) > params.eccLim(1))  ) );
    sigmaXFail = (params.applySigma== true ...
	&&(    (currCluster.ddparams(sColX) < params.sigmaLim(1)) ||(currCluster.ddparams(sColX) > params.sigmaLim(2)) ...
	    || (currCluster.aaparams(sColX) < params.sigmaLim(1)) ||(currCluster.aaparams(sColX) > params.sigmaLim(2)) ...
	    || (currCluster.daparams(sColX) < params.sigmaLim(1)) ||(currCluster.daparams(sColX) > params.sigmaLim(2)) ) );
    sigmaYFail = (params.applySigma== true ...
	&&(    (currCluster.ddparams(sColY) < params.sigmaLim(1)) ||(currCluster.ddparams(sColY) > params.sigmaLim(2)) ...
	    || (currCluster.aaparams(sColY) < params.sigmaLim(1)) ||(currCluster.aaparams(sColY) > params.sigmaLim(2)) ...
	    || (currCluster.daparams(sColY) < params.sigmaLim(1)) ||(currCluster.daparams(sColY) > params.sigmaLim(2)) ) );
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
DDlinkedpos=zeros(0,2);
AAlinkedpos=zeros(0,2);
DAlinkedpos=zeros(0,2);

nClusterOut = numel(filteredCluster);
for i = 1:nClusterOut
  DDlinkedpos=[DDlinkedpos;filteredCluster(i).dd];
  AAlinkedpos=[AAlinkedpos;filteredCluster(i).aa];
  DAlinkedpos=[DAlinkedpos;filteredCluster(i).da];
end



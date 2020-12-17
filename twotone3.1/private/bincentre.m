function [xcentre,xedge] = bincentre(xminedge,xmaxedge,nbin)
%function [xcentre,xedge] = bincentre(xminedge,xmaxedge,nbin)
% calculate bin centres for n equally spaced bins defined by upper and lower edges xminedge and xmaxedge
% 
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 Alpha, released 101115
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%
nedge = nbin + 1;
xedge = linspace(xminedge,xmaxedge,nedge);
xcentre = xedge(1:end-1)+diff(xedge)/2;


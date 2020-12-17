function detectedParticles = findpeaks3(imfilt,im, threshold,  bpDiameter,windowSize)
%function detectedParticles = findpeaks3(imfilt,im, threshold,  bpDiameter,windowSize)
% plot the points using the emory routines
% modified 241008 to take account of bandpass setting edges 
%
% outputs a matrix containing:
%	| xPos | yPos | brightness |sx | sy| ecc |
%
%  SH 100418 modified to remove the size exclusion threshold
% 
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 Alpha, released 101115
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%

posLim =0.5*windowSize;%if its going to move significantly we already have a problem
sigmaLim = [0 50];
%make sure the image is a double
im = double(im);

% Locate the peaks to pixel level accuracy
out =findLocalMax(imfilt, threshold);
if isempty(out)
	detectedParticles = zeros(0,6);
else
	if numel(out) >0
		% Locate the centroids to sub pixel level accuracy on the original unfilterd image
		%tic
		[points brightness sx sy ecc]= lsqLoc(im,out(:,1:2),windowSize,posLim,sigmaLim);
		%toc
	else
		points = out;
	end

	detectedParticles = zeros(size(points,1),6);

	detectedParticles = [points(:,1:2), brightness, sx, sy, ecc];
end
%ECC = ecc;
%DETECTEDPARTICLES = detectedParticles;
%figure; hist(ecc,[0.01:0.02:0.99]);
%----------------------------------------------------------------------
function out=findLocalMax(im,th)
% finds local maxima in an image to pixel level accuracy.   
%  local maxima condition is >= rather than >
% inspired by Crocker & Griers algorithm, and also Daniel Blair & Eric Dufresne's implementation
%   im = input image for detection
%   th - detection threshold - pixels must be strictly greater than this value
% out : x,y coordinates of local maxima

%identify above threshold pixels
[y,x] = find(im>th);
%delete pixels identified at the boundary of the image
[imsizey, imsizex]= size(im);
edgepixel_idx = find( y==1 | y==imsizey | x==1 | x==imsizex);
y(edgepixel_idx) = [];
x(edgepixel_idx) = [];

%check if its a local maxima
subim = zeros(3,3);
islocalmaxima = zeros(numel(x),1);
for i = 1:numel(x)
  subim = im([y(i)-1:y(i)+1],[x(i)-1:x(i)+1]);
  islocalmaxima(i) = all(subim(2,2)>=subim(:));
end
%assign all above threshold pixels to out initially
out = [x,y];
%delete pixels which are not local maxima
out(find(~islocalmaxima),:) = [];
%--------------------------------------------
function [posOutAll brightness sx sy ecc] = lsqLoc(im,pos,windowSize,posLim,sigmaLim);
%function [posOutAll brightness sx sy ecc] = lsqLoc(im,pos,windowSize,posLim,sigmaLim);

if exist('isGaussFitToolsInstalled','file') && isGaussFitToolsInstalled()==1
  useCPPfit = true;
  fitArg = {};
else
  useCPPfit = false;
  warning('GaussFitTools PSF fitting library not installed - expect extremely slow performance!');
  fitArg = {'useMatlabFit'};
end

nPos = size(pos,1);
k = 1;
for i =1:nPos
  posCur = pos(i,:);
  initguess = [0];%no initguess
  %initguess = [amplitude widthguess background X_POSim Y_POSim ];

  [phot_count a normChi2 posOut eccentricity] = freeGaussFitEllipse( im, pos(i,:), windowSize,posLim, sigmaLim, initguess,'autoDetect',fitArg{:});

  if ~all(a == 0 )%if the fit has worked
    posOutAll(k,:) = posOut;
    brightness(k,1) = a(1);
    sx(k,1) = a(2);
    sy(k,1) = a(3);
    ecc(k,1) = eccentricity;
    k=k+1;
  end
end

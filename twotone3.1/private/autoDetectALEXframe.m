function firstGreenFrame = autoDetectALEXframe(moviePath,imageLim, varargin)
%function firstGreenFrame = autoDetectALEXframe(moviePath,imageLim, varargin)
% function firstGreenFrame = autoDetectALEXframe(moviePath,imageLim)
% function firstGreenFrame = autoDetectALEXframe(moviePath,imageLim, ..., 'framesToAnalyse', nFrames, ...)
% function firstGreenFrame = autoDetectALEXframe(moviePath,imageLim, ..., 'skipFirstFrame', ...)
%
% Autodetect whether the first green frame is the first or 2nd frame in the movie
% OPTIONAL PARAMETERS:
%	 'framesToAnalyse' parameter is how many frames to calculate the average based on - default is 10
%	 'skipFirstFrame' skips the first frame in the calculation - since this is often blank for typical - default is true
%	 data on our setups
% 
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 Alpha, released 101115
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only�? license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%

% get the even & odd frame averages
[oddGreenFrameAvg, evenGreenFrameAvg, oddRedFrameAvg, evenRedFrameAvg] = ... 
		getFrameAverages(moviePath, imageLim, varargin{:});

% calculate the totals
oddGreenFrameSum = sum(oddGreenFrameAvg(:));
evenGreenFrameSum = sum(evenGreenFrameAvg(:));
%oddRedFrameSum = sum(oddRedFrameAvg(:));
%evenRedFrameSum = sum(evenRedFrameAvg(:));

%if (oddGreenFrameSum > evenGreenFrameSum) && (evenRedFrameSum > oddRedFrameSum)
%	firstGreenFrame = 1;
%elseif  (oddGreenFrameSum < evenGreenFrameSum) && (evenRedFrameSum < oddRedFrameSum)
%	firstGreenFrame = 2;
%else
%	display('red & green averages give conflicting information');
%	error('autoDetectALEXframe:cannotCalculateALEXframe' ,'Unable to unambiguously determine first green frame');
%end
disp(oddGreenFrameSum)
disp(evenGreenFrameSum)
if (oddGreenFrameSum > evenGreenFrameSum) 
	firstGreenFrame = 1;
elseif  (oddGreenFrameSum < evenGreenFrameSum)
	firstGreenFrame = 2;
end
% -------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------
function [oddGreenFrameAvg, evenGreenFrameAvg, oddRedFrameAvg, evenRedFrameAvg] = getFrameAverages(moviePath,imageLim, varargin);
% function [oddGreenFrameAvg, evenGreenFrameAvg, oddRedFrameAvg, evenRedFrameAvg] = getFrameAverages(moviePath,imageLim, varargin);
%
% get the even & odd frame averages

% set defaults
% number of frames to analyse
nFramesAnalyse= 8;
% skipFirstFrame
skipFirstFrame = true;

%parse optional arguments
m = nargin;
n = length(varargin);
i = 1;
while i <= n
	if strcmp(varargin{i}, 'skipFirstFrame')
		skipFirstFrame = varargin{i+1};
		i = i + 2;
	elseif strcmp(varargin{i}, 'framesToAnalyse')
		nFramesAnalyse= varargin{i+1};
		if isnan(nFramesAnalyse)... % is not a number
		   || rem(nFramesAnalyse, 1)~=0 ... % is not integer
		   || (nFramesAnalyse < 1)
			error('nFramesAnalyse must be integer >= 1');
		end
		i = i + 2 ;
	else
		error('unrecognised option supplied');
		i = i + 1;
	end
end

tirfIm = TirfImage(moviePath,0,imageLim);
greenStack = getGreenStack(tirfIm);
redStack = getRedStack(tirfIm);

numFrames = getNumFrames(greenStack);

%set the frames to analyse over
if skipFirstFrame == false
	startFrame = 1;
	firstOddFrame = 1;
elseif numFrames == 2
	startFrame = 1;
	firstOddFrame = 1;
else
	startFrame = 2;
	firstOddFrame = 3;
end

% basic error checking
if startFrame == numFrames
	error('Movie is too short for alex analysis - only 1 frame contains data');
end

lastFrame = startFrame + nFramesAnalyse - 1;
if lastFrame > numFrames 
	lastFrame = numFrames;
end

% initialise the avg matrices
oddGreenFrameAvg = zeros(size(getFrame(greenStack, 1)));
evenGreenFrameAvg = oddGreenFrameAvg;
oddRedFrameAvg = zeros(size(getFrame(redStack, 1)));
evenRedFrameAvg = oddRedFrameAvg;

% get the frame averages
numOddFramesAvg = numel(firstOddFrame:2:lastFrame);
numEvenFramesAvg = numel(2:2:lastFrame);

for i = firstOddFrame:2:lastFrame
	greenFrame = double( getFrame(greenStack, i));
	oddGreenFrameAvg =  oddGreenFrameAvg + greenFrame/numOddFramesAvg;

	redFrame = double( getFrame(redStack, i));
	oddRedFrameAvg =  oddRedFrameAvg + redFrame/numOddFramesAvg;
end

for i = 2:2:lastFrame
	greenFrame = double( getFrame(greenStack, i));
	evenGreenFrameAvg =  evenGreenFrameAvg + greenFrame/numEvenFramesAvg;

	redFrame = double( getFrame(redStack, i));
	evenRedFrameAvg =  evenRedFrameAvg + redFrame/numEvenFramesAvg;
end


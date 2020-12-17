function ImStack = ImageStack(impath, imageLim )
% ImageStack class constructor
% 
%	ImStack = ImageStack(impath, imageLim )
% 
% Description: 
%	Creates an object allowing access to the channels
% of a given tirfim frame %
% Inputs:
%	impath
%	imageLim is a 1x4 matrix containing the limits of the channel
%	  ie, [ Xstart, Xend, Ystart, Yend];
%	
% Data:	
%	ImPath
%	ImageLim
%	NumFrames
%	ImType - this is the numeric data type of the image eg 'uint16'
%	FileType - this is the file type - either 'fits' or 'tif'
%	NAXIS1 :size x of whole movie
%	NAXIS2 :size y ofwhole movie
%	ImageInfo - stack containing information about the image
% Public Methods:
%	- get[all]
%	- getFrame( ImStack, n)
% Private Methods:
%
% History:
% 	130208 	- First alpha complete (SH)
%	180208 	- Function headers and help standardised (SH)
%	200208  - Completely revised to make it more memory efficient
%	210708  - Takes the numFramesAvg as an opt. argument to only average the first n frames
%
% Author:
%	Seamus Holden
%
%
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 Alpha, released 101115
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only�? license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%
%

% work out the file type
[pathstr, name, ext] = fileparts(impath);
if strcmp(ext,'.fits')
  ImStack.FileType = 'fits';
elseif strcmp(ext,'.tif') || strcmp(ext,'.tiff')
  ImStack.FileType = 'tif';
else
  error('Unrecognised file type.');
end

if strcmp(ImStack.FileType,'fits')
  %ImageInfo = fits_read_header(impath);
  ImageInfo  = fitsinfo(impath);
  keywords = ImageInfo.PrimaryData.Keywords; 
  IndexC = strfind(keywords(:,1), 'NAXIS3');
  Index = find(not(cellfun('isempty', IndexC)));
  %check whether there is more than one frame
  if ~isempty(Index)
	  ImStack.NumFrames= keywords{Index,2};
  else
	  ImStack.NumFrames= 1;
  end

  %these are the widths and the heights need them to read fits frames
  IndexC = strfind(keywords(:,1), 'NAXIS1');
  Index = find(not(cellfun('isempty', IndexC)));
  ImStack.NAXIS1 = keywords{Index,2};
  IndexC = strfind(keywords(:,1), 'NAXIS2');
  Index = find(not(cellfun('isempty', IndexC)));
  ImStack.NAXIS2 = keywords{Index,2};
  ImStack.ImageInfo = ImageInfo;
  ImStack.movie = fitsread(impath);
elseif strcmp(ImStack.FileType,'tif')
  ImageInfo = imfinfo(impath);
  ImStack.NumFrames=numel(ImageInfo);
  ImStack.NAXIS1 = ImageInfo(1).Width;
  ImStack.NAXIS2 = ImageInfo(1).Height;
  ImStack.ImageInfo = ImageInfo;
else
  error('Unrecognised file type.');
end

%check the supplied image limits are within bounds
if (imageLim(1) < 1) || (imageLim(2) < 1) ...
      || (imageLim(1) > ImStack.NAXIS1) || (imageLim(2) > ImStack.NAXIS1) ...
      || (imageLim(1) >= imageLim(2) ) ...
      ||(imageLim(3) < 1) || (imageLim(4) < 1) ...
      || (imageLim(3) > ImStack.NAXIS2) || (imageLim(4) > ImStack.NAXIS2) ...
      || (imageLim(3) >= imageLim(4) )
  warning('ImageStack:OutOfBounds','Supplied image limits exceed the bounds of the image, using whole image!');
  imageLim(1) = 1;imageLim(2)=ImStack.NAXIS1;
  imageLim(3) = 1;imageLim(4)=ImStack.NAXIS2;
end 

ImStack.ImageLim = imageLim;

ImStack.ImPath = impath;

ImStack.ImType = NaN;
% Setup the class
ImStack = class(ImStack, 'ImageStack');

ImStack.ImType = class(getFrame(ImStack,1));
%ImType - this is the numeric data type of the image eg 'uint16' form an average image
% need to convert to double temporarily to avoid saturation 
% and then convert back to the correct format


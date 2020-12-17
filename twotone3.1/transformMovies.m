function transformMovies()
% transformMovies()
% 
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 , released 110426
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%
% FUNCTION: transformMovies
% DESCRIPTION:
%   carry out image transforms and beads calibration
%
%   You manually pick the control points using matlabs cpselect function
%   The donor channel is presented on the left, acceptor channel on the right.
%   Once you are finished, click "File">"Close Control Point Selection Tool"
%   and you will be presented with a transformed overlay of the two images
%   If they look correct, accept the transform and it will be saved.
%
%   We use a projective transform (read cp2tform documentation)
%   If you want to change the transform type, you can modify the "align" function
%   starting at line ~246.
% Returned channel is the transform which maps donor molecules into the acceptor channel
% ie 
%     x_A  = T x_D (in matrix notation
% You apply this transform by running 
%   x_A = tformfwd(TFORM,x_D);
% in matlab
% to go the other way, just run
%   x_D = tformfwd(TFORM,x_A);
%
%

FIRSTGREENFRAME = 0 ; %CW image always supplied for transform
ALTERNATIONPERIOD = 2; %Assume were dealing with twocolour
AVGFIRST = 1;
AVGLAST  = 20;

display(' ');
display('----TirfImage calibration and transformation----')

if ~exist('beadspath','var')
  % Get the calibration image
  input('Press return to load calibration image file');
  
  filterindex = 0;
  while filterindex == 0 %check user has supplied a filename
      [filename,pathname, filterindex] = uigetfile({'*.fits;*.tif'; '*.*'});
  end
  
  beadspath = strcat( pathname,filename);
end

%load the current image limits
if exist(twotoneInstallDir('','twotoneDefaultSettings_CW.mat'),'file')
  load(twotoneInstallDir('','twotoneDefaultSettings_CW.mat'),'twotoneData');
  imageLim = twotoneData.settings.imageSettings.channelImageLim;
else %create twotone settings file
  initializeTwotoneSetting('CWonly')
  load(twotoneInstallDir('','twotoneDefaultSettings_CW.mat'),'twotoneData');
  imageLim = twotoneData.settings.imageSettings.channelImageLim;
end

%check whether we need to update the image limits
fprintf('\nCurrent image limits are:\n');
fprintf('Donor:\t\t');
fprintf('%d ',imageLim(1,:));
fprintf('\n');
fprintf('Acceptor:\t');
fprintf('%d ',imageLim(2,:));
fprintf('\n');
changelim = getAnswer('Would you like to change the image limits (y/N) : ', 'y', 'n','','Y','N');
if strcmp(changelim, 'y')||strcmp(changelim, 'Y')
  imageLimOK = 'n';
  while strcmp(imageLimOK,'n')

    imageLim = zeros(2,4);
    fprintf('Please input pixel limits for each emmision channel\n')
    imageLim(1,1) = getIntegerAnswer('Donor X start: ' , 1, inf);
    imageLim(1,2) = getIntegerAnswer('Donor X end: ' , 1, inf);
    imageLim(1,3) = getIntegerAnswer('Donor Y start: ' , 1, inf);
    imageLim(1,4) = getIntegerAnswer('Donor Y end: ' , 1, inf);
    imageLim(2,1) = getIntegerAnswer('Acceptor X start: ' , 1, inf);
    imageLim(2,2) = getIntegerAnswer('Acceptor X end: ' , 1, inf);
    imageLim(2,3) = getIntegerAnswer('Acceptor Y start: ' , 1, inf);
    imageLim(2,4) = getIntegerAnswer('Acceptor Y end: ' , 1, inf);
    fprintf('Pixel limits inputs were: \n');
    fprintf('Donor:\t\t');
    fprintf('%d ',imageLim(1,:));
    fprintf('\n');
    fprintf('Acceptor:\t');
    fprintf('%d ',imageLim(2,:));
    fprintf('\n');

    imageLimOK= getAnswer('Is this correct (y/n) : ', 'y', 'n');
  end
  
  twotoneData.settings.imageSettings.channelImageLim = imageLim;
  save(twotoneInstallDir('','twotoneDefaultSettings_CW.mat'), 'twotoneData');
  fprintf('Updated image limits saved sucessfully.\n');
end

fprintf('\nLoading file %s\n',beadspath);
% Load the image as a CalImage object
calIm = TirfImage( beadspath, FIRSTGREENFRAME, imageLim, ALTERNATIONPERIOD);

% Get the control points from the user
modify_cal = 'y';
firstLoop = 1;
while strcmp(modify_cal,'y')
  
  % Initiate the user calibration
  if firstLoop == 1
    [TFORM, green_in_points,red_base_points]= calibrate(calIm, AVGFIRST, AVGLAST);
    firstLoop = 0;
  else
    [TFORM, green_in_points,red_base_points]= calibrate(calIm, AVGFIRST, AVGLAST,green_in_points,red_base_points);
  end

  modify_cal = getAnswer('Do you want to modify the calibration? (y/n)', 'y', 'n');

end
if ~exist('TFORMpath','var')  
  % save the TFORM file
  display(' ');
  display('Please select a path to save the transform file:');
  
  filterindex = 0;
  while filterindex == 0 %check user has supplied a filename
    [filename,pathname, filterindex] = uiputfile('*.tform.mat', ...
              'Save TFORM');
  end
  
  TFORMpath = strcat(pathname, filename);
  
  %If the file doesnt have the correct ending, append it
  if ~strcmp(TFORMpath(end-9:end),'.tform.mat') 
    
    if strcmp(TFORMpath(end-3:end),'.mat')
      TFORMpath = strcat(TFORMpath(1:end-4),'.tform.mat');
    else
      TFORMpath = strcat(TFORMpath,'.tform.mat');
    end
  
  end
end

save( TFORMpath, 'TFORM');
display('TFORM saved to file'); 
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
function answer = getAnswer(question, varargin)
%   answer = getAnswer(Question, varargin)
% function to stop having to write loops to get a quick 
% user input

num_options = nargin - 1;

sensible_answer = false;

% loop round the input options to check if youve got an answer
while sensible_answer == false
  choice = input(question, 's');

  for i = 1:num_options
  
    if strcmp(choice, varargin{i}) == 1
      sensible_answer = true;
      answer = choice;
    end
  end
  
  if sensible_answer == false
    display('Not a valid choice');
  end
end
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
function answer = getIntegerAnswer(question, lowerbound, upperbound)
%function answer = getIntegerAnswer(question, lowerbound, upperbound)
% function to stop having to write loops to get a quick 
% integer user input

sensible_answer = false;

% loop round the input options to check if youve got an answer
while sensible_answer == false
  choice = input(question, 's');
  choice = str2double(choice);

  if isnumeric(choice) && rem(choice , 1) == 0 ...
      && choice>=lowerbound && choice<=upperbound
    sensible_answer = true;
    answer = choice;
  end
  
  if sensible_answer == false
    display('Not a valid choice');
  end
end
%---------------------------------------------------------------
function [TFORM, green_in_points,red_base_points] =calibrate( calIm,AVGFIRST, AVGLAST,green_in_points,red_base_points)
%function [TFORM, green_in_points,red_base_points] =calibrate( calIm,AVGFIRST, AVGLAST,green_in_points,red_base_points)
% obtains an image transform by manual alignment of user, from 
% an input calibration image - of the CalImage class  

firstGreenFrame = getFirstGreen(calIm);

[DexDem,AexAem, DexAem, framesData] = calculateAveragesALEX(calIm, firstGreenFrame, AVGFIRST, AVGLAST);

%transformtype = 'projective';
% get the inputs for the alignment function
green_in = DexDem;
red_base = AexAem;

if ~exist('green_in_points','var')
  green_in_points = NaN;
  red_base_points = NaN;
end

% align the green to the red
% outTFORM is your transform, *_points are your control points
% default assignment of control points is NaN, ie before the first calibration
% This is to flag up an empty control point set.
% if statement check whether any control points already exist

% The try statement is error catching to ensure enough points are picked
% Instead of quitting, the tool loops until enough points are picked

repeat_align = 1;

while repeat_align ==1
  try 
    if any(any(isnan(green_in_points))) || any(any(isnan(red_base_points)))
      [green_registered, TFORM, green_in_points, red_base_points] = align(green_in, red_base);
    else
      [green_registered, TFORM, green_in_points, red_base_points] = align(green_in, red_base, green_in_points, red_base_points);
    end
    repeat_align= 0;
  catch ME
    if  ME.identifier(1:23)=='Images:cp2tform:atLeast'
      warning('Not enough points to infer transform.');
      repeat_align = 1;
    else
      rethrow(ME);
    end
  end
end

% Update the CalIm control points
green_points = green_in_points;
red_points = red_base_points;

%-----------------------------------------------------------------------------------------------
function [green_registered, mytform, green_input_points, red_base_points] = align(green_input, red_base, green_input_points, red_base_points)
%[green_registered, mytform,green_input_points, red_base_points] = 
%  align(green_input, red_base, green_input_points, red_base_pointsmsho)
% Produce a transform matrix for TIRF FRET analysis 
% beads alignment image
% ouputs green image adjusted to red & transform matrix
% given red & green beads files

% adjust the dynamic range so that the whole of the space 0->65535 is used
% see stretchlim help for detaile
% TODO Allow different types of transform130808
%for i = 1:nargin
%  if 
%transformtype = 'polynomial';
transformtype = 'projective';
%polynomial order
%polOrder = 4;

%normalise the images
green_input = double(green_input);
green_input = (green_input - min(green_input(:))) ./ ( max(green_input(:)) - min(green_input(:))) ;
red_base = double(red_base);
red_base = (red_base - min(red_base(:))) ./ ( max(red_base(:)) - min(red_base(:))) ;

[sizey sizex] = size(red_base);
%set the control points using the cp select tool
%if previous cps are supplied use these
keepLooping = true;
while keepLooping == true
  if nargin == 2
    [green_input_points, red_base_points] = cpselect(green_input, red_base, 'Wait', true);
  elseif nargin ==4
    [green_input_points, red_base_points] = cpselect(green_input, red_base, green_input_points, red_base_points, 'Wait', true);
  else
    error('incorrect number of arguments (2 or 4)');
  end
  if isempty(green_input_points)
    warning('At least 4 points required to generate transform');
    keepLooping = true;
  else	
    keepLooping = false;
  end
end

%calculate global transform using the selected points
if exist('polOrder')
  mytform = cp2tform(green_input_points, red_base_points, transformtype, polOrder);
else
  mytform = cp2tform(green_input_points, red_base_points, transformtype);
end

%apply global transform to green image
green_registered = imtransform(green_input, mytform, 'FillValues', 0, 'xdata', [1 sizex], 'ydata', [1 sizey]);

% overlay the red and green_registered images for visual inspection
%by creating a dummy true colour image of the data
dummy_blue = zeros(size(red_base));
visual = cat(3, red_base, green_registered, dummy_blue);
figure; h1 = imshow(visual);

clear visual dummy_blue;


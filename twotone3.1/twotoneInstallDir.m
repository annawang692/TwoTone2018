function pathInsideInstallDir = twotoneInstallDir(varargin)
% pathInInstallDir = twotoneInstallDir(varargin)
% 
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 , released 110426
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%
% FUNCTION: twotoneInstallDir
% return the twotone install path, along with any required
% appended directories and files supplied in varargin
%	Example:
%  pathInsideInstallDir = twotoneInstallDir('defaultConfigFiles', '1x1binning.autoDetConfig.mat')
% 	returns
%  pathInsideInstallDir =

% path of currentFunction
currentFunctionPath = mfilename('fullpath');
% this is in the install directory 
twotoneInstallDirectory = fileparts(currentFunctionPath);
% movie must be in the same directory as the positions file
if nargin == 1 % ie just a directory
	pathInsideInstallDir = fullfile(twotoneInstallDirectory, varargin{:});
	pathInsideInstallDir = strcat(pathInsideInstallDir, filesep);
elseif nargin == 2 % ie a filename is supplied too

	pathInsideInstallDir = fullfile(twotoneInstallDirectory, varargin{:});
else
	pathInsideInstallDir = fullfile(twotoneInstallDirectory, '');
	pathInsideInstallDir = strcat(pathInsideInstallDir, filesep);
end



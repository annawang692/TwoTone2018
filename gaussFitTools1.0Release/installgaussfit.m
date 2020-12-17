function installtwotone()
% linux users who have matlab installed to a shared location will need to
% run this as sudo.

ostype = computer();
% path of currentFunction
currentFunctionPath = mfilename('fullpath');
% this is in the install directory 
twotoneBundleDirectory = fileparts(currentFunctionPath);

gaussLib = 'GaussFitTools';

if strcmp(ostype,'PCWIN')
  gaussBin = 'GaussFitTools_bin_Win32';
elseif strcmp(ostype,'GLNX86')%linux 32bit 
  gaussBin = 'GaussFitTools_bin_Linux32';
elseif strcmp(ostype,'GLNXA64')%linux 64bit 
  gaussBin = 'GaussFitTools_bin_Linux64';
end

gausspath   = fullfile(twotoneBundleDirectory,gaussLib,gaussBin);

addpath(gausspath);
savepath();


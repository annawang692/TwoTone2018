function [logText] = twotoneAnalyzeCW(fileName,twotoneData);
% function [logText] = twotoneAnalyzeCW(fileName,twotoneData);
% 
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 Alpha, released 101115
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only�? license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%

logText = [];

% get the name of the movie
[path name extension] = fileparts(fileName);

%set some info about the data
twotoneData.results.movieInfo.path  = path;
twotoneData.results.movieInfo.name  = fileName;
twotoneData.results.analysisDate    = datestr(now,'yymmdd');
% work out the file type
if strcmp(extension,'.fits')
  ImageInfo = fitsread(fileName);
elseif strcmp(extension,'.tif') || strcmp(ext,'.tiff')
  ImageInfo = imfinfo(fileName);
else
  error('Unrecognised file type.');
end
twotoneData.results.movieInfo.imageInfo  = ImageInfo;
  
%1. carry out autodetection
autoDetectIsOK = false;
try
  logText= sprintfappend('For file %s, generating auto-detected positions file\n', fileName);
  fprintf('For file %s, generating auto-detected positions file\n', fileName);
  
  % run autodetection
  [twotoneData, clusteredData, logTextAutoDet] = autoDetectCW(fileName,twotoneData);
  logText = [logText,logTextAutoDet];

  autoDetectIsOK = true;

catch ME %if firstGreenFrame auto detection fails, make a note of it 
    %in the log file but dont crash
    
  if strcmp(ME.identifier, 'autoDetectALEXframe:cannotCalculateALEXframe')
    %auto-detect alternation failed - make a note & move onto next file
    %print to screen
    fprintf('Auto-detection of the alternation order for ALEX movie failed\n');
    fprintf('for file %s\n', fileName);
    fprintf('Skipping this file\n\n');
    %print to file
    logText= sprintfappend(logText,'Auto-detection of the alternation order for ALEX movie failed\n');
    logText= sprintfappend(logText,'for file %s\n', fileName);
    logText= sprintfappend(logText,'Skipping this file\n\n');
  elseif strcmp(ME.identifier,'autoDetectMain:TFORMerror')
         %Transform file does not exist
    %print to screen
    fprintf('Transform file does not exist for file %s\n', fileName);
    fprintf('Skipping this file\n\n');
    %print to file
    logText= sprintfappend(logText,'Transform file does not exist for file %s\n', fileName);
    logText= sprintfappend(logText,'Skipping this file\n\n');
  else
    %print to screen
    fprintf('Autodetection failed for unknown reason for file %s. Error report is:\n\n', fileName);
    fprintf('%s',getReport(ME));
    fprintf('\n');
    fprintf('Skipping this file\n\n');
    %print to file
    logText= sprintfappend(logText,'Autodetection failed for unknown reason for file %s. Error report is:\n\n', fileName);
    logText= sprintfappend(logText,'%s',getReport(ME));
    logText= sprintfappend(logText,'\n');
    logText= sprintfappend(logText,'Skipping this file\n\n');
  end
end

%do fitting
if autoDetectIsOK == true
  try
    %fit the data
    fprintf('For file %s, extracting data\n', fileName);
    logText= sprintfappend(logText,'For file %s, extracting data\n', fileName);
    [twotoneData] = extractDataCW(fileName,twotoneData,clusteredData);
    %logText = [logText,logTextExtractData];
    %save the data
    savePath = [fileName, twotoneData.settings.fileAppendString];
    fprintf('Saving data to %s\n\n', savePath);
    logText= sprintfappend(logText,'Saving data to %s\n\n', savePath);
    save(savePath,'twotoneData');
  catch ME
    if strcmp(ME.identifier,'MATLAB:minrhs')
      fprintf('No points to analyse for file\n');
      logText= sprintfappend(logText, 'No points to analyse for file\n');
      fprintf('for file %s\n', fileName);
      logText= sprintfappend(logText, 'for file %s\n', fileName);
      fprintf('Skipping this file\n\n');
      logText= sprintfappend(logText, 'Skipping this file\n\n');
    else
      fprintf('Data extraction failed for unknown reason for file %s. Error report is:\n\n', fileName);
      logText= sprintfappend(logText, 'Data extraction failed for unknown reason for file %s. Error report is:\n\n', fileName);
      fprintf('%s',getReport(ME));
      logText= sprintfappend(logText, '%s',getReport(ME));
      fprintf('\n');
      logText= sprintfappend(logText, '\n');
      fprintf('Skipping this file\n\n');
      logText= sprintfappend(logText, 'Skipping this file\n\n');
    end
  end
end



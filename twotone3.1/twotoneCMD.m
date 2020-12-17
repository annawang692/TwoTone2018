function logFilePath = twotoneCMD(fileFilter, configPath, varargin)
% logFilePath = twotoneCMD(fileFilter, configPath)
% 
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0, released 110426
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only�? license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%
% FUNCTION: twotoneCMD
% DESCRIPTION: Analyse TIRF-FRET data in batch mode.%
% INPUTS:
%   - fileFilter:  string containg an expression for which all matching files in 
%     local directory will be analysed
%  - configPath: path to the config file used to carry out the analysis. 
% Optional arguments:
%  'OutputTextData', true/false 
%  'OutputTwotoneOldFormat', true/false 
% OUTPUTS:
%  - logFilePath: path of the log file summarising results of analysis
%

%parse the optional arguments
%  'OutputTextData', true/false 
%  'OutputTwotoneOldFormat', true/false 
outputTextData = false;
outputTwotoneOldFormat = false;
n = length(varargin);
i = 1;
while i <= n 
  if strcmp(varargin{i},'OutputTextData')
    outputTextData = varargin{i+1};
    i = i + 2;
  elseif strcmp(varargin{i},'OutputTwotoneOldFormat')
    outputTwotoneOldFormat = varargin{i+1};
    i = i + 2;
  else
    error('unexpected argument supplied');
  end
end


load(configPath,'twotoneData');

%set the path of the log file
CURRENTTIMESTRING = datestr(clock, '.HHMM.ddmmyy');
logFilePath = strcat('batchAnalysisLog',CURRENTTIMESTRING,'.txt');

%Initialise the log file
fid = fopen(logFilePath, 'w+');
%record the setup data
fprintf(fid,'Time: %s\n', datestr(clock,'HHMM ddmmyy') );
fprintf(fid, 'Beginning batch analysis of TIRF images\n\n');
fprintf(fid,'Using config file:\n');
fprintf(fid,'%s\n', configPath);
fprintf(fid,'Config file contents:\n');

%write same info to screen
fprintf('Time: %s\n', datestr(clock,'HHMM ddmmyy') );
fprintf('Beginning batch analysis of TIRF images\n');
fprintf('Using config file:\n');
fprintf('%s\n', configPath);

fprintf(fid,'Beginning analysis for files matching\n');
if iscell(fileFilter)
  fprintf(fid,'''%s''\n', fileFilter{:});
else
  fprintf(fid,'''%s''\n', fileFilter);
end

if  twotoneData.settings.imageSettings.alternationPeriod == 2
  fprintf(fid,'Beginning ALEX analysis\n');
  fprintf('Beginning ALEX analysis\n');
elseif twotoneData.settings.imageSettings.alternationPeriod == 1
  fprintf(fid,'Beginning CW analysis\n');
  fprintf('Beginning CW analysis\n');
 else
    error('Unrecognised alternation period');
end

%parse filter
[files, num_files] = parseFileFilter(fileFilter);
%filesOut counter, n
n = 1;
filesOut = {};
for i = 1:num_files
  
  currentFilePath = files{i};
  if twotoneData.settings.imageSettings.alternationPeriod == 2
    [logText] = twotoneAnalyzeALEX(files{i},twotoneData);
  elseif  twotoneData.settings.imageSettings.alternationPeriod == 1
    [logText] = twotoneAnalyzeCW(files{i},twotoneData);
  else
    error('Unrecognised alternation period');
  end
  logText = strrep(logText,'\','\\');
  fprintf(fid, logText);
  if outputTextData == true
    try
      analyzedFileName = [files{i},twotoneData.settings.fileAppendString];
      convertTwotoneToText(analyzedFileName,[analyzedFileName,'.txt']);
    catch ME
      fprintf('\n Conversion to text output failed.\n');
    end
  end

  if outputTwotoneOldFormat == true
    try
      analyzedFileName = [files{i},twotoneData.settings.fileAppendString];
      convertTwotoneToOldFormat(analyzedFileName,[analyzedFileName,'.v2.mat']);
    catch ME
      fprintf('\n Conversion to v2 output format failed.\n');
    end
  end
end

%close the file
fclose (fid);

fid = fopen(logFilePath, 'a');
fprintf(fid, '\n\n%%-------------------------------------------------------------------\n');
fprintf(fid, '%%-------------------------------------------------------------------\n\n');
fprintf(fid,'\n\nBatch analysis completed sucessfully\n');
%record the time
currentTime = datestr(clock,'HHMM ddmmyy') ;
fprintf(fid,'Time: %s\n', currentTime );

fclose(fid);

% print the same info to the screen
fprintf( '\n\n%%-------------------------------------------------------------------\n');
fprintf( '%%-------------------------------------------------------------------\n\n');
fprintf('\n\nBatch analysis completed sucessfully\n');
%record the time
fprintf('Time: %s\n', currentTime );
%-------------------------------------------------------------------
%-------------------------------------------------------------------
%-------------------------------------------------------------------
end

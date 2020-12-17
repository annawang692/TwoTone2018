function [files, numFiles] = parseFileFilter(fileFilter);
% function [files, numFiles] = parseFileFilter(fileFilter);
%
% take the input to the batch processing functions, 'fileFilter'
% work out whether its a string (a regex fileFilter matching several files
% or whether its a cell containing filenames
% return filenames which match it and number of files
% 
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 Alpha, released 101115
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%

if iscell(fileFilter) % if its a cell of file names
	numFiles = numel(fileFilter);
	files = fileFilter;
elseif isstr(fileFilter) % if its a string
	
	filesDummy = dir(fileFilter);
	numFiles = numel(filesDummy);
	files = {};
	for i = 1:numFiles
		files{i} = filesDummy(i).name;
	end

else
	error('could not parse supplied fileFilter');
end


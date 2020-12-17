%compile script
%
% MFITSIO2 Version 1.0, author S Holden, University of Oxford
% DERIVED FROM MFITSIO 1.2.3, author Damian Eads, LNL.
% For licensing information, see 'COPYING'
% For file authorship details, see ChangeLog
%

file ={...
'fits_delete_keyword.c'
'fits_read_header.c'
'fits_read_image.c'
'fits_read_image_subset.c'
'fits_write_header.c'
'fits_write_image.c'
'fits_write_image_subset.c'
};

ostype = computer();

if strcmp(ostype,'PCWIN')%win 32bit install
  binFolder = 'mfitsio2_bin_Win32';
  cflags = {'-Dsnprintf=_snprintf',['-L',binFolder],['-I',binFolder],'-lcfitsio'};%the linux library snprintf is called _snprintf in windows
  c_files = {'mfitsio.c'};
  for i = 1:numel(file)
    mex(file{i},c_files{:},cflags{:});
  end
  outfileType = 'mexw32';
  movefile(['*',outfileType,'*'],binFolder);
  copyfile(['*.m'],binFolder);
elseif strcmp(ostype,'GLNXA64')%linux 64bit install - requires install of libcfitsio3 & libcfitsio3-dev from repositories
  binFolder = 'mfitsio2_bin_Linux64';
  cflags = {'-lcfitsio','-I/usr/local/include'};
  c_files = {'mfitsio.c'};
  for i = 1:numel(file)
    mex(file{i},c_files{:},cflags{:});
  end
  outfileType = 'mexa64';
  movefile(['*',outfileType,'*'],binFolder);
  copyfile(['*.m'],binFolder);
elseif strcmp(ostype,'GLNX86')%linux 32bit install - requires install of libcfitsio3 & libcfitsio3-dev from repositories
  binFolder = 'mfitsio2_bin_Linux32';
  cflags = {'-lcfitsio','-I/usr/local/include'};
  c_files = {'mfitsio.c'};
  for i = 1:numel(file)
    mex(file{i},c_files{:},cflags{:});
  end
  outfileType = 'mexglx';
  movefile(['*',outfileType,'*'],binFolder);
  copyfile(['*.m'],binFolder);
end


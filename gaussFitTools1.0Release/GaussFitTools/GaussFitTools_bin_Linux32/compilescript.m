%compile script
%
%GaussFitTools 
%Copyright (C) 2010 The Chancellor, Masters and Scholars of the University of Oxford.
%Authors: Oliver J Britton, Seamus J Holden.
%
%This program is free software; you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation; either version 2 of the License, or
%(at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%

file ={...
'gaussfit_fixedposition.cpp',
'gaussfit_fixedposition_elliptical.cpp',
'gaussfit_fixedsigma.cpp',
'gaussfit_fixedsigmaposition.cpp',
'gaussfit_free.cpp',
'gaussfit_free_elliptical.cpp',
};

ostype = computer();

if strcmp(ostype,'PCWIN')%win 32bit install
  binFolder = 'GaussFitTools_bin_Win32';
  for i = 1:numel(file)
    mex(file{i},'-L', ['-L',binFolder], '-lblas_win32', '-llapack_win32', '-llevmar')
  end
  outfileType = 'mexw32';
  movefile(['*',outfileType,'*'],binFolder);
  copyfile(['*.m'],binFolder);
elseif strcmp(ostype,'GLNXA64')%linux 64bit install 
  FLAGS_64BIT='-fPIC';
  binFolder='GaussFitTools_bin_Linux64';

  for i = 1:numel(file)
    mex('-O',['CFLAGS=',FLAGS_64BIT],['-L',binFolder], '-llevmar','-lgfortran',file{i}, [binFolder,'/lapack_LINUX.a'], [binFolder,'/blas_LINUX.a'])
  end
  outfileType = 'mexa64';
  movefile(['*',outfileType,'*'],binFolder);
  copyfile(['*.m'],binFolder);
elseif strcmp(ostype,'GLNX86')%linux 32 bit install - requires liblapack-dev libblas-dev
  binFolder='GaussFitTools_bin_Linux32';
  c_flags = {['-L',binFolder],  '-llevmar','-llapack', '-lblas' }; 
  for i = 1:numel(file)
    mex(file{i}, c_flags{:})
  end
  outfileType = 'mexglx';
  movefile(['*',outfileType,'*'],binFolder);
  copyfile(['*.m'],binFolder);

end


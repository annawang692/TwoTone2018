% Function Name:
%    fits_write
%
% Description: Writes an image and header to a FITS file. 'BITPIX', 'NAXIS',
% and 'SIMPLE' keywords are ignored. This information is automatically
% calculated based on the dimensions and the data type of the input image.
%
% Usage:
%    fits_write(FILENAME, HEADER, IMAGE);
%
% Arguments:
%    FILENAME: A character array representing the filename.
%    HEADER:   The fits header as a structure array.
%    IMAGE:    The image as a MATLAB array.
%
% Type 'mfitsio_license' to display the MFITSIO licensing agreement.

function fits_write(FILENAME, HEADER, IMAGE)
  fits_write_image(FILENAME, IMAGE, HEADER);
%  fits_write_header(FILENAME, HEADER);

% MFITSIO2 Version 1.0, author S Holden, University of Oxford
% DERIVED FROM MFITSIO 1.2.3, author Damian Eads, Los Alamos National Laboratory.
% For licensing information, see 'COPYING'
% For file authorship details, see ChangeLog

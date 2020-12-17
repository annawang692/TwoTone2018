% Function Name:
%    fits_read
%
% Description: Reads a FITS file and stores all header information in
% a structure array.
%
% Usage:
%    [HEADER, IMAGE] = fits_read(FILENAME);
%
% Arguments:
%    FILENAME: A character array representing the filename.
%
% Returns:
%    HEADER:   The fits header as a structure array.
%    IMAGE:    The image as a MATLAB array.
%
% Type 'mfitsio_license' to display the MFITSIO licensing agreement.

function [HEADER, IMAGE] = fits_read(filename);
HEADER = fits_read_header(filename);
IMAGE = fits_read_image(filename);

% MFITSIO2 Version 1.0, author S Holden, University of Oxford
% DERIVED FROM MFITSIO 1.2.3, author Damian Eads, LNL.
% For licensing information, see 'COPYING'
% For file authorship details, see ChangeLog


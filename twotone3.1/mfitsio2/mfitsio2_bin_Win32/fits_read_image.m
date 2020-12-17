% Function Name:
%    fits_read_image
%
% Description: Reads an image from a FITS file.
%
% Usage:
%    [IMAGE] = fits_read_image(FILENAME);
%
% Arguments:
%    FILENAME: A character array representing the filename.
%
% Returns:
%    IMAGE:    The image as a MATLAB array.
%
% Type 'mfitsio_license' to display the MFITSIO licensing agreement.

function [IMAGE]=fits_read_image(FILENAME);

% MFITSIO2 Version 1.0, author S Holden, University of Oxford
% DERIVED FROM MFITSIO 1.2.3, author Damian Eads, LNL.
% For licensing information, see 'COPYING'
% For file authorship details, see ChangeLog

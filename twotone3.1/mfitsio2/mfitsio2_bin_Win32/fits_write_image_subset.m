% Function Name:
%    fits_write_image_subset
%
% Description: Writes an image to a region of a FITS file image.
%
% Usage:
%    [IMAGE] = fits_write_image_subset(FILENAME, IMG, START);
%    [IMAGE] = fits_write_image_subset(FILENAME, IMG, START, HEADER);
%
% Arguments:
%    FILENAME: A character array representing the filename.
%    IMG:      The image to write to the region of the FITS image.
%    START:    A vector representing the starting position of the region 
%              to start writing the input image.
%    HEADER:   Optional. A structure array containing keywords to be modified.
%
% Returns:
%    IMAGE:    The image as a MATLAB array.
%
% Type 'mfitsio_license' to display the MFITSIO licensing agreement.

function fits_write_image_subset(FILENAME, IMG, START, HEADER);

% MFITSIO2 Version 1.0, author S Holden, University of Oxford
% DERIVED FROM MFITSIO 1.2.3, author Damian Eads, LNL.
% For licensing information, see 'COPYING'
% For file authorship details, see ChangeLog


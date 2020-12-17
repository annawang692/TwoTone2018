% Function Name:
%    fits_read_image_subset
%
% Description: Reads a region of an image from a FITS file. This function
%              is particularly useful for programs which must process
%              large images.
%
% Usage:
%    [IMAGE] = fits_read_image_subset(FILENAME, START, END);
%
% Arguments:
%    FILENAME: A character array representing the filename.
%    START:    A vector representing the starting position of the image region.
%    END:      A vector representing the ending position of the image region.
%
% Returns:
%    IMAGE:    The image as a MATLAB array.
%
% Type 'mfitsio_license' to display the MFITSIO licensing agreement.

function [IMAGE]=fits_read_image_subset(FILENAME, START, END);

% MFITSIO2 Version 1.0, author S Holden, University of Oxford
% DERIVED FROM MFITSIO 1.2.3, author Damian Eads, LNL.
% For licensing information, see 'COPYING'
% For file authorship details, see ChangeLog


% Function Name:
%    fits_write_image
%
% Description: Writes an image and a header to a FITS file.
%
% Usage:
%    fits_write_image(FILENAME, IMAGE);
%    fits_write_image(FILENAME, IMAGE, HEADER);
%
% Arguments:
%    FILENAME: A character array representing the filename.
%    IMAGE: An image array representing the filename.
%    HEADER: Optional. A structure array containing keywords to be modified.
%
% Returns:
%    Nothing.
%
% Type 'mfitsio_license' to display the MFITSIO licensing agreement.

function fits_write_image(FILENAME, IMAGE, HEADER);

% MFITSIO2 Version 1.0, author S Holden, University of Oxford
% DERIVED FROM MFITSIO 1.2.3, author Damian Eads, Los Alamos National Laboratory.
% For licensing information, see 'COPYING'
% For file authorship details, see ChangeLog


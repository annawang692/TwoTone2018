% Function Name:
%    fits_write_header
%
% Description: Write a header to a FITS file. Any 'BITPIX', 'NAXIS', or
% 'SIMPLE' field is ignored. The function does not rewrite the entire header.
% Instead, it only rewrites fields present in the HEADER structure.
%
% Usage:
%    fits_write_header(FILENAME, HEADER);
%
% Arguments:
%    FILENAME: A character array representing the filename.
%    HEADER: A structure array containing the header keywords to modify.
%
% Returns:
%    Nothing.
%
% Type 'mfitsio_license' to display the MFITSIO licensing agreement.

function fits_write_header(FILENAME, HEADER);

% MFITSIO2 Version 1.0, author S Holden, University of Oxford
% DERIVED FROM MFITSIO 1.2.3, author Damian Eads, Los Alamos National Laboratory.
% For licensing information, see 'COPYING'
% For file authorship details, see ChangeLog

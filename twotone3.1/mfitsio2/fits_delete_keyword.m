% Function Name:
%    fits_delete_keyword
%
% Description: Deletes a keyword from the header of a FITS file.
%
% Usage:
%    fits_delete_keyword(FILENAME, KEYWORD);
%
% Arguments:
%    FILENAME: A character array representing the filename.
%    KEYWORD: The keyword to delete.
%
% Returns:
%    Nothing.
%
% Type 'mfitsio_license' to display the MFITSIO licensing agreement.

function fits_delete_keyword(FILENAME, KEYWORD);
% MFITSIO2 Version 1.0, author S Holden, University of Oxford
% DERIVED FROM MFITSIO 1.2.3, author Damian Eads, Los Alamos National Laboratory.
% For licensing information, see 'COPYING'
% For file authorship details, see ChangeLog

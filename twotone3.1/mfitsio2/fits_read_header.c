/**
 * File:        fits_read_header.c
 * Author:      Damian Ryan Eads <eads@lanl.gov>
 * Description: The gateway program for fits_read_header.
 * Created:     December 2, 2002
 *
% MFITSIO2 Version 1.0, author S Holden, University of Oxford
% DERIVED FROM MFITSIO 1.2.3, author Damian Eads, LNL.
% For licensing information, see 'COPYING'
% For file authorship details, see ChangeLog
 */

#include <string.h>
#include <stdio.h>
#include <fitsio.h>
#include <malloc/malloc.h>

#include "mex.h"
#include "matrix.h"
#include "mfitsio.h"

/**
 * Executes the fits_read_header function. Type "help fits_read_header" for
 * more information.
 *
 * @param nlhs             The total number of output arguments.
 * @param plhs             The output arguments.
 * @param nrhs             The total number of input arguments.
 * @param prhs             The input arguments.
 */

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  int buflen, status;
  char *filename;
  mfitsio_header *header;

  /* Check proper input and output */
  if (nrhs != 1)
    mexErrMsgTxt ("One input required.");
  else if (nlhs > 1)
    mexErrMsgTxt ("Too many output arguments.");
  else if (!mxIsChar (prhs[0]))
    mexErrMsgTxt ("Input must be a string (filename).");

  if (mxGetM (prhs[0]) != 1)
    mexErrMsgTxt ("Input must be a row vector string (filename).");

  buflen = (mxGetM (prhs[0]) * mxGetN (prhs[0])) + 1;

  filename = (char*) mxCalloc (buflen, sizeof (char));
  status = mxGetString (prhs[0], filename, buflen);

  header = mfitsio_read_header (filename);
  if (header != 0)
  {
    plhs[0] = mfitsio_adapt_fheader (header);
    mfitsio_free_header (header);
  }
  else
  {
    mexErrMsgTxt ("Unable to open FITS file.");
  }

  return;
}

/**
 * File:        fits_read_image.c
 * Author:      Damian Ryan Eads <eads@lanl.gov>
 * Description: The gateway program for fits_read_image.
 * Created:     December 2, 2002
 *
% MFITSIO2 Version 1.0, author S Holden, University of Oxford
% DERIVED FROM MFITSIO 1.2.3, author Damian Eads, Los Alamos National Laboratory.
% For licensing information, see 'COPYING'
% For file authorship details, see ChangeLog
 */

#include <string.h>
#include <stdio.h>
#include <fitsio.h>
#include </System/Library/Frameworks/Kernel.framework/Versions/A/Headers/sys/malloc.h>

#include "mex.h"
#include "matrix.h"
#include "mfitsio.h"

/**
 * Executes the fits_read_image function. Type "help fits_read_image"
 * for more information.
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
  mfitsio_info *info;
  /**  mxArray *headerm;**/

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

  filename = (char*)mxCalloc (buflen, sizeof (char));
  status = mxGetString (prhs[0], filename, buflen);

  header = mfitsio_read_header (filename);
  if (header != 0)
  {
    /**    headerm = mfitsio_adapt_fheader (header);**/
    info = mfitsio_read_info (filename);
    plhs[0] = mfitsio_read_image (filename, info);
    mfitsio_free_header (header);
    mfitsio_free_info (info);
  }
  else
  {
    mexErrMsgTxt ("Unable to open FITS file.");
  }

  return;
}

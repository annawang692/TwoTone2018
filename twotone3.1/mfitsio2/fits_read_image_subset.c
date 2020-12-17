/**
 * File:        fits_read_image_subset.c
 * Author:      Damian Ryan Eads <eads@lanl.gov>
 * Description: The gateway program for fits_read_image_subset.
 * Created:     December 2, 2002
 * Changes: 	S Holden 2009
% MFITSIO2 Version 1.0, author S Holden, University of Oxford
% DERIVED FROM MFITSIO 1.2.3, author Damian Eads, LNL.
% For licensing information, see 'COPYING'
% For file authorship details, see ChangeLog
 * Modified by SH 090204 to fix a memory leak
 * The leak is somewhere in the accessing of the header file - to correct it
 * remove all accesses of the header file
 */


#include <string.h>
#include <stdio.h>
#include <fitsio.h>
#include <malloc/malloc.h>

#include "mex.h"
#include "matrix.h"
#include "mfitsio.h"

/**
 * Executes the fits_read_image_subset function. Type
 * "help fits_read_image_subset" for more information.
 *
 * @param nlhs             The total number of output arguments.
 * @param plhs             The output arguments.
 * @param nrhs             The total number of input arguments.
 * @param prhs             The input arguments.
 */

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  int buflen, status, i;
  char *filename;
  long *fpixels, *lpixels, *lnaxes, *inc;
  mfitsio_info *info;
  const mxArray *startp, *endp;
  double *dataStart, *dataEnd;

  /* Check proper input and output */
  if (nrhs != 3)
    mexErrMsgTxt ("Three inputs required.");
  else if (nlhs > 1)
    mexErrMsgTxt ("Too many output arguments.");
  else if (!mxIsChar (prhs[0]))
    mexErrMsgTxt ("Input must be a string (filename).");

  if (mxGetM (prhs[0]) != 1)
    mexErrMsgTxt ("Input must be a row vector string (filename).");

  buflen = (mxGetM (prhs[0]) * mxGetN (prhs[0])) + 1;

  filename = (char*) mxCalloc (buflen, sizeof (char));
  status = mxGetString (prhs[0], filename, buflen);
	
  /*change here - use info to check if the file exists instead of header */
  info = mfitsio_read_info (filename);

  startp = prhs[1];
  endp = prhs[2];

  /** If the header exists. */
  if (info != 0)
  {

    /** Retrieve a pointer to the data of the starting and ending pixels.*/

    mfitsio_check_coordinate(startp, endp, info->naxis, info->naxes);
    dataStart = (double*) mxGetData(startp);
    dataEnd = (double*) mxGetData(endp);
    fpixels = mfitsio_convert_dbl2long(dataStart, info->naxis);
    lpixels = mfitsio_convert_dbl2long(dataEnd, info->naxis);
    inc = mfitsio_create_ones_vector(info->naxis);
    lnaxes = mfitsio_convert_int2long(info->naxes, info->naxis);
    plhs[0] = mfitsio_read_image_impl (filename, info, fpixels, lpixels,
				       lnaxes, inc);
    mfitsio_free_info (info);
    MFITSIO_FREE(fpixels);
    MFITSIO_FREE(lpixels);
    MFITSIO_FREE(inc);
    MFITSIO_FREE(lnaxes);
  }
  else
  {
    mexErrMsgTxt ("Unable to open FITS file.");
  }

  return;
}

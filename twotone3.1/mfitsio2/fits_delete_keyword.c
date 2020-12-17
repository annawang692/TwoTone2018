/**
 * File:        fits_delete_keyword.c
 * Author:      Damian Ryan Eads <eads@lanl.gov>
 * Description: The gateway program for fits_delete_keyword.
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
#include <malloc/malloc.h>

#include "mex.h"
#include "matrix.h"
#include "mfitsio.h"

/**
 * Executes the fits_delete_keyword function. Type "help fits_delete_keyword"
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
  char *filename, *keyword;
  mfitsio_header *header;

  /* Check proper input and output */
  if (nrhs != 2)
  {
    mexErrMsgTxt ("Two inputs required.");
  }
  else if (nlhs > 1)
  {
    mexErrMsgTxt ("Too many output arguments.");
  }
  else if (!mxIsChar (prhs[0]) || !mxIsChar (prhs[1]))
  {
    mexErrMsgTxt ("Inputs must be strings (filename, keyword).");
  }
  else if (mxGetM (prhs[0]) != 1 || mxGetM (prhs[1]) != 1)
  {
    mexErrMsgTxt ("Inputs must be a row vector strings (filename, keyword).");
  }
  buflen = (mxGetM (prhs[0]) * mxGetN (prhs[0])) + 1;
  filename = (char *) mxCalloc (buflen, sizeof (char));
  status = mxGetString (prhs[0], filename, buflen);

  buflen = (mxGetM (prhs[1]) * mxGetN (prhs[1])) + 1;
  keyword = (char *) mxCalloc (buflen, sizeof (char));
  status = mxGetString (prhs[1], keyword, buflen);


  mfitsio_delete_keyword (filename, keyword);
  return;
}

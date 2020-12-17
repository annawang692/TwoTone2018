/**
 * File:        fits_write_image.c
 * Author:      Damian Ryan Eads <eads@lanl.gov>
 * Description: The gateway program for fits_write_image.
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
 * Executes the fits_write_image function. Type "help fits_write_image" for
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
  mfitsio_info *info;
  mxArray *headerm = 0;
  mxArray *img = 0;

  /* Check proper input and output */
  if (nrhs < 2)
    mexErrMsgTxt ("At least two inputs required.");
  else if (nlhs > 0)
    mexErrMsgTxt ("Too many output arguments.");
  else if (!mxIsChar (prhs[0]))
    mexErrMsgTxt ("Input must be a string (filename).");

  if (mxGetM (prhs[0]) != 1)
    mexErrMsgTxt ("Input must be a row vector string (filename).");

  buflen = (mxGetM (prhs[0]) * mxGetN (prhs[0])) + 1;

  filename = (char*) mxCalloc (buflen, sizeof (char));
  status = mxGetString (prhs[0], filename, buflen);
  img = (mxArray *) prhs[1];

  /**  if (mxGetClassID(img) != mxDOUBLE_CLASS) {
    mexErrMsgTxt("The image should be a M x N x P x ... double array. Use "
                 "double(IMG) to cast.");
		 }**/
  if (nrhs == 3)
  {
    headerm = (mxArray *) prhs[2];
    if (mxGetClassID (headerm) != mxSTRUCT_CLASS)
    {
      mexErrMsgTxt ("The header must be a 1x1 struct array.");
    }
  }
  else
  {
    int dims[] = { 1, 1 };
    int ndims = 2;
    int nfields = 1;
    const char **fields;
    fields = (const char **) malloc (sizeof (char *) * nfields);
    fields[0] = "_nil_";
    headerm = mxCreateStructArray (ndims, dims, nfields, fields);
  }
  mfitsio_write_image (filename, headerm, img);

  return;
}


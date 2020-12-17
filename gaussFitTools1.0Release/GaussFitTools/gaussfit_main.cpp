#include "mex.h"

void gaussfit_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], char* variant, char* model, const int m, const int background_int,
	void(*fit_model)(double *p, double *x, int m, int n, void *data), void(*fit_jacobian)(double *p, double *jac, int m, int n, void *data),
	const int numEqualityConstraints, const int numBounds, const int numInputs, const int numDataEntries)

{

	//gets the dimensions of the input data (window size)
	int sizex = mxGetN(prhs[0]);
	int sizey = mxGetM(prhs[0]);


	//checks that input frame is square
	if (sizex != sizey)
	{
		mexErrMsgTxt("Input data is not a square matrix.");
	}
	//===================================


	//Read in data from MATLAB
	//===================================

	int n = sizex*sizey;    //number of pixels in data
							//here we create a dynamic array of size N to hold the data

	double* x;
	x = new double[n];


	//Read in data to x from MATLAB
	double *datapoint = mxGetPr(prhs[0]);

	for (int j = 0; j < sizey; j++)
	{
		for (int i = 0; i < sizex; i++)
		{

			x[(i * sizex) + j] = datapoint[i + (j*sizey)];
		}
	}

	//Because of the way MATLAB reads mxArrays (down each column)
	//important to account for this because of having fixed positions etc
	//or fixed elliptical width.

	//so read the MATLAB array consecutively for speed and then the first element (1,1) goes to
	//position (1,1) [element 0], then the second element (1,2) goes to position (2,1) [element 1*sizex] etc
	//so we fill x[] down the "columns" in the array
	//eg for a 3x3 matrix we read in elements[0 - 8] from MATLAB and input them in array element order
	//[0], [3], [6] (first column done), [1],[4],[7] (second), [2], [5], [8] - all done

	double* constraints_pointer;

	double *data;           //read in data from MATLAB input, if required
	if (numDataEntries == 0)
	{
		data = NULL;
	}
	else
	{
		data = new double[numDataEntries];

		for (int i = 0; i < numDataEntries; i++)
		{
			constraints_pointer = mxGetPr(prhs[i + 1]); //+1 because array of intensity data is indexed at [0]
			data[i] = constraints_pointer[0];
		}

	}
	//data contains any constants we need to pass to the constraints functions
	//data contains any constants we need to pass to the constraints functions


	//Read in bounding condition data from MATLAB to pass to constraints function

	//Matlab input convention - data first, then all fixed parameters, then all bounded parameters in order -
	//lower bound first, then upper bound, with standard variable sequence A, sigma, b, Xo, Yo, then extra options. So bound indices for
	// prhs[] start at numEqualityConstraints + 1 ([0] is always the data)

	//===================================


	//Sets up initial guess and boundaries for solver
	//===================================
	double* p;
	p = new double[m];
	double back_p[1];

	void initialfit(char *variant, double p[], double data_array[], int n, double *data); //define initial fit function

																						  //if flag in input is set to true, use supplied guess, otherwise run initial fit routine
	bool* use_initial_fit_guess = mxGetLogicals(prhs[numInputs - 3]);
	if (use_initial_fit_guess[0] == true) //then use the guess given in input
	{
		double* init_guess_pointer = mxGetPr(prhs[numEqualityConstraints + numBounds + 2]); //gets pointer to guess
		for (int i = 0; i < m; i++)
		{
			p[i] = init_guess_pointer[i];
		} //transfer guess from input argument


	}

	else //run the initial fits program to generate the initial guess
	{
		initialfit(variant, p, x, n, data);
	}

	//sets the initial value for the background to be used by the background only solver
	back_p[0] = p[background_int - 1];

	//initialise boundary and penalty weighting arrays
	//which will be taken as the arguments for constraint calculating
	//and solver functions - constraints() and dlevmar_etc()

	double* lb;
	double* ub;
	lb = new double[m];
	ub = new double[m];
	//lower and upper boundary arrays

	double* weights;
	weights = new double[m]; //Not actually used but if we want to define weighting for penalty function
							 //for use in the boundary conditions, this is here


	double info[LM_INFO_SZ];

	//Read in search stopping parameters from input array from Matlab

	double* options_pointer = mxGetPr(prhs[numEqualityConstraints + numBounds + 3]);

	double searchstopParams[4] = { options_pointer[1], options_pointer[2], options_pointer[3] };
	int numIterations = options_pointer[4];

	//this array controls parameters epsilon 1,2 and 3 as described in http://www.ics.forth.gr/~lourakis/levmar/levmar.pdf
	//the purpose of them is to stop the solver before it has reached its maximum number of iterations if the solution has
	//already sufficiently converged
	//for data with reasonable signal to noise, sufficient convergence should only take ~10 iterations, so these conditions
	//are important

	//searchstopParams[0] is epsilon1, the gradient of chi squared, 10^-2 or -3 seems to work well
	//searchstopParams[1] is epsilon2, the gradient of the change in best fit parameters, relative to the values of those parameters
	//this is generally the condition that will determine the speed/accuracy/precision of the output and is worth playing around with
	//higher values generally slow down the routine, but increase the accuracy, especially for badly scaled conditions (eg. background
	//several orders of magnitude higher than peak amplitude) good values seem to be in the range 10^-3 (for well scaled conditions) to
	//10^-5
	//searchstopParams[2] is epsilon3, the absolute value of chi squared. For a well solved gaussian with low noise
	//this will still be a large number so this is essentially ignored





	double opts[LM_OPTS_SZ] = { LM_INIT_MU, searchstopParams[0],searchstopParams[1],searchstopParams[2], LM_DIFF_DELTA };
	//opts[0] is initial mu value - controls the damping in the LM algorithm - using default
	//opts[1] is stopping threshold on gradient of chi square
	//opts[2] is stopping threshold on relative change of parameter magnitudes
	//opts[3] is stopping threshold on value of chi squared
	//opts[4] is irrelevant if analytic jacobian used, so use default
	//===================================

	//define constraint setting function
	void constraints(char *variant, double *lb, double *ub, double *constraintMatrix, double *constraintRHS,
		int numBounds, int numEqualityConstraints, double *weights, int m, int framesize, double *data);

	int ret;








	//We will always have box (inequality constraints) eg. amplitude must be positive
	//However, in some cases we will have linear (equality) constraints and will need to use dlevmar_blec_der
	//But, if we have no equality constraints we should use dlevmar_bc_der which uses box constraints only
	//So use an if-else condition on numEqualityConstraints to decide which to us

	if (numEqualityConstraints > 0) //then we need box and linear constraints

	{
		double *ConstraintMatrix = new double[numEqualityConstraints * m];
		double *ConstraintRHS = new double[numEqualityConstraints];

		///constraints function
		constraints(variant, lb, ub, ConstraintMatrix, ConstraintRHS, numBounds, numEqualityConstraints, weights, m, sizex, data);

		///main solver function if you need equality constraints
		ret = dlevmar_blec_der

			///basics
			(fit_model,  //pointer to the function that calculates the model values at each datapoint given parameters p
				fit_jacobian, //pointer to the function that calculates the analytic jacobian dxi/dpj at each datapoint
				p, //current best guess parameters as an array of size m
				x,  //data to be fitted to the model, array of size n
				m,  //number of free parameters
				n, //number of datapoints

				   ///constraints
				lb, //array of m elements setting lower bounds on each parameter (p[i] >= lb[i])
				ub, //array of m elements setting upper bounds on each parameter (p[i] <= ub[i])
				ConstraintMatrix, //matrix of equality constraint values, numEqualityconstraints = number of rows
								  //m is number of columns
				ConstraintRHS, //RHS of constraint equation ConstraintMatrix * p = ConstraintRHS
				numEqualityConstraints, //how many constraints you have
				NULL, //penalty weighting array, set to NULL to use defaults, otherwise fill and use weights[m]

					  ///solver stop procedure and debugging
				numIterations, //max no. of iterations the solver will run for, 50 to 100 seems about right
				opts, //5 element array, described above
				info, //from source code of levmar library:
					  /* info: information regarding the minimization. Used to find the lsq minimisation value at output
					  * info[0]= ||e||_2 at initial p.
					  * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||Dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
					  * info[5]= # iterations,
					  * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
					  *                                 2 - stopped by small Dp
					  *                                 3 - stopped by itmax
					  *                                 4 - singular matrix. Restart from current p with increased mu
					  *                                 5 - no further error reduction is possible. Restart with increased mu
					  *                                 6 - stopped by small ||e||_2
					  *                                 7 - stopped by invalid (i.e. NaN or Inf) "func" values. This is a user error
					  * info[7]= # function evaluations
					  * info[8]= # Jacobian evaluations
					  * info[9]= # linear systems solved, i.e. # attempts for reducing error
					  */

				NULL,    // working memory at least LM_BLEC_DER_WORKSZ() reals large, allocated if NULL
				NULL,   //Covariance matrix corresponding to LS solution; mxm. Set to NULL if not needed.
				data); // pointer to additional data passed to function and jacobian calculator - fixed constants etc

		delete[] ConstraintMatrix;
		delete[] ConstraintRHS;
	}

	///main solver function if you only need inequality constraints - see above for arg details
	else if (numEqualityConstraints == 0) //we only need box constraints

	{
		///constraints function, inequality constraints only
		constraints(variant, lb, ub, NULL, NULL, numBounds, numEqualityConstraints, weights, m, sizex, data);
		///main solver, see above for arg descriptions
		ret = dlevmar_bc_der

		(fit_model, fit_jacobian, p, x, m, n,

			lb, ub,
			
			NULL, numIterations, opts, info, NULL, NULL, data); // with analytic Jacobian



	}

	else { mexErrMsgTxt("Number of equality constraints not set properly, should be an integer >= 0"); }





	///OUTPUT
	//output brightness, final fit parameters and normalised least squares value to matlab

	double brightness, lsq_normalised;
	const double pi = 3.14159265358979323846; //easier to do this than require a library Pi function

	bool got_signal = true; //assume signal as default

	if (model == "circular")
	{
		brightness = 2.0*pi*p[0] * pow(p[1], 2);
	}
	else if (model == "elliptical")
	{
		brightness = 2.0*pi*p[0] * p[1] * p[2];
	}
	else
	{
		mexErrMsgTxt("Error, model not set.");
	}

	//if elliptical model transform the paramaters st. s_y is semimajor axis, and 0<theta<pi
	if (model == "elliptical")
	{
		double s_x, s_y, theta;
		s_x = p[1];
		s_y = p[2];
		theta = p[6];

		if (s_x > s_y) //ie if the fit has got the axes the wrong way around
		{
			p[1] = s_y;//swap the two values
			p[2] = s_x;
			theta = theta + pi / 2; //rotate the axes 90deg
		}

		//transform theta st. its 0<theta<pi ie 
		//ie, mod(theta,pi) = theta - floor(theta/pi) * pi;

		theta = theta - floor(theta / pi) * pi;
		p[6] = theta;
	}




	lsq_normalised = info[1] / n;



	plhs[0] = mxCreateDoubleScalar(brightness); //send brightness as output to MATLAB

	plhs[1] = mxCreateDoubleMatrix(1, m, mxREAL);

	//need to do this next bit because can only pass an array back to Matlab as a matrix if array is
	//in dynamic memory via mxCalloc command
	void* final_fit_alloc = mxCalloc(m, sizeof(double));
	double* final_fit = (double *)final_fit_alloc;


	for (int i = 0; i < m; i++)
	{
		final_fit[i] = p[i];
	}

	mxFree(mxGetPr(plhs[1])); //mxSetPr will not deallocate memory allocated to plhs[1], so do it manually or get a leak
	mxSetPr(plhs[1], final_fit); //fit parameters are now ready to be passed to Matlab


	plhs[2] = mxCreateDoubleScalar(lsq_normalised);


	if (options_pointer[0] == true) //then display the output, otherwise just forward it to MATLAB
	{
		if (model == "circular")
		{
			mexPrintf("\nBrightness is %f.\n", brightness);
			mexPrintf("Fit parameters are: A = %f, sigma = %f, b = %f, Xo = %f, Yo = %f. \n", p[0], p[1], p[2], p[3], p[4]);
			mexPrintf("Normalised lsq is %f.\n", lsq_normalised);
		}

		if (model == "elliptical") //then display the output, otherwise just forward it to MATLAB
		{
			mexPrintf("\nBrightness is %f.\n", brightness);
			mexPrintf("Fit parameters are: A = %f, sigma = %f, %f, b = %f, Xo = %f, Yo = %f. Theta = %f. \n", p[0], p[1], p[2], p[3], p[4], p[5], p[6]);
			mexPrintf("Normalised lsq is %f.\n", lsq_normalised);
		}
		//******************Uncomment for extra information*************************************** 
		mexPrintf("Levenberg-Marquardt returned in %g iter, reason %g, sumsq %g [%g]\n", info[5], info[6], info[1], info[0]);
		//   Debugging - info[0] returns lsq at start of fit,info[1] returns lsq at end of fit,
		//******************Uncomment for extra information*************************************** 
	}

	//freeup memory allocated for dynamic arrays
	delete[] x;
	delete[] p;
	delete[] data;
	delete[] lb;
	delete[] ub;
	delete[] weights;


	return;


}

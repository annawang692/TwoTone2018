//Constraint setting functions
/***%
%GaussFitTools 
%Copyright (C) 2010 The Chancellor, Masters and Scholars of the University of Oxford.
%Authors: Oliver J Britton, Seamus J Holden.
%
%This program is free software; you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation; either version 2 of the License, or
%(at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
*/

void constraints(char *variant, double *lb, double *ub, double *constraintMatrix, double *constraintRHS, int numBounds, int numEqualityConstraints, double *weights, int m, int framesize, double *data)


//function sets up constraints
//lb and ub do less than or equal to and greater than or equal to constraints on each of the m parameters

//constraintMatrix (M) and Rhs set up a system of constraint equations M p = rhs

//weights sets some kind of penalty function thing - not understood currently


//FOR WRITING NEW CONSTRAINTS - put new function definition below this comment,
//run them in the else if section below, and add the function to the bottom of this file
{
    void
    constraints_basic(double *lb, double *ub, int numBounds, int numEqualityConstraints, double *weights, int m, int framesize, double* data);

    void
    constraints_fixedposition(double *constraintMatrix, double *constraintRHS, int numBounds, int numEqualityConstraints, int m, double* data);

    void
    constraints_fixedsigma(double *constraintMatrix, double *constraintRHS, int numBounds, int numEqualityConstraints, int m, double* data);

    void
    constraints_fixedsigmaposition(double *constraintMatrix, double *constraintRHS, int numBounds, int numEqualityConstraints, int m, double* data);

    void
    constraints_basic_elliptical(double *lb, double *ub, int numBounds, int numEqualityConstraints, double *weights, int m, int framesize, double *data);

    void
    constraints_fixedposition_elliptical(double *constraintMatrix, double *constraintRHS, int numBounds, int numEqualityConstraints, int m, double* data);
    //declare new functions here


    //PROGRAM BEGINS



    if(variant == "free")
        {
            constraints_basic(lb,ub,numBounds,numEqualityConstraints,weights,m,framesize,data);
            //Sets the standard inequality constraints that you will (probably) always want:
            //Positive A, non-negative sigma smaller than framesize, positive b, psf centered within frame
            //and constraints specified by matlab inpu

            return;
        }

    else if(variant == "fixedposition")
        {
            constraints_basic(lb,ub,numBounds,numEqualityConstraints,weights,m,framesize,data);
            constraints_fixedposition(constraintMatrix, constraintRHS, numBounds, numEqualityConstraints, m, data);
            return;
        }

    else if(variant == "fixedsigma")
        {
            constraints_basic(lb,ub,numBounds,numEqualityConstraints,weights,m,framesize,data);
            constraints_fixedsigma(constraintMatrix, constraintRHS, numBounds, numEqualityConstraints, m, data);
            return;
        }

    else if(variant == "fixedsigmaposition")
        {
            constraints_basic(lb,ub,numBounds,numEqualityConstraints,weights,m,framesize,data);
            constraints_fixedsigmaposition(constraintMatrix, constraintRHS, numBounds, numEqualityConstraints, m, data);
            return;
        }

    else if(variant == "free_elliptical")
        {
            constraints_basic_elliptical(lb,ub,numBounds,numEqualityConstraints,weights,m,framesize,data);
            return;
        }

    else if(variant == "fixedposition_elliptical")
        {
            constraints_basic_elliptical(lb,ub,numBounds,numEqualityConstraints,weights,m,framesize,data);
            constraints_fixedposition_elliptical(constraintMatrix, constraintRHS, numBounds, numEqualityConstraints, m, data);
            return;
        }

    else
        {
            mexErrMsgTxt("Variant error - unknown variant specified please check char* variant = "" at top of gaussmexfits_...");
        }


  //include new functions here


        return;
}


//basic bounding constraints for a circular gaussian, bounds for sigma and position specified by the user
void constraints_basic(double *lb, double *ub, int numBounds, int numEqualityConstraints, double *weights, int m, int framesize, double *data)
{


/// bounds - A,b hardcoded, others specified from Matlab input

//A constraints  - A is greater than or equal to zero (non-negative)
lb[0] = 0.0;
ub[0] = 65535.0;


//sigma constraints
lb[1] = data[numEqualityConstraints]; //equality constraints are loaded into data first, so the bounds come after
ub[1] = data[numEqualityConstraints+1];

//b constraints - b takes any value
lb[2] = -65535.0;
ub[2] = 65535.0;

//Xo and Yo constraints
lb[3] = data[numEqualityConstraints+2];
ub[3] = data[numEqualityConstraints+3]; //Xo

lb[4] = data[numEqualityConstraints+4];
ub[4] = data[numEqualityConstraints+5]; //Yo



}

///FIXED Xo, Yo (position)
void constraints_fixedposition(double *constraintMatrix, double *constraintRHS, int numBounds, int numEqualityConstraints, int m, double* data)

    {

        for (int i = 0; i < numEqualityConstraints; i++)  //fill constraintMatrix and RHS row by row
         {

             for (int j = 0; j < m; j++)
                {

                 constraintMatrix[j + (i * m)] = 0.0;

                }
            constraintRHS[i] = 0.0;

         } //fill each array with zeros

   double Xo = data[0];
   double Yo = data[1];
                               //(first 0->numBounds-1 of data contain the bounding constraints
    //Now want constraints Xo = data[numBounds] Yo = data[numBounds+1]. constraintMatrix is 2x5
    //So only interested in elements (1,4) and (2,5) - 4th and 10 elements (levmar uses row major order)
    constraintMatrix[3] = 1;
    constraintRHS[0] = Xo; //Xo = input

    constraintMatrix[9] = 1;
    constraintRHS[1] = Yo; //Yo = input

    }

///FIXED SIGMA
void constraints_fixedsigma(double *constraintMatrix, double *constraintRHS, int numBounds, int numEqualityConstraints, int m, double* data)

    {
        for (int i = 0; i < numEqualityConstraints; i++) //fill constraintMatrix and RHS row by row
         {

             for (int j = 0; j < m; j++)
                {

                 constraintMatrix[j + (i * m)] = 0.0;

                }
            constraintRHS[i] = 0.0;

         } //fill each array with zeros

   double sigma = data[0];

    //Now want constraint sigma = data[0], constraintMatrix is 1x5
    //So only interested in element (1,2), the 2nd element(levmar uses row major order)
    constraintMatrix[1] = 1;
    constraintRHS[0] = sigma; //sigma = input

    }



///FIXED SIGMA AND POSITION
void constraints_fixedsigmaposition(double *constraintMatrix, double *constraintRHS, int numBounds, int numEqualityConstraints, int m, double* data)

    {
        for (int i = 0; i < numEqualityConstraints; i++)  //fill constraintMatrix and RHS row by row
         {

             for (int j = 0; j < m; j++)
                {

                 constraintMatrix[j + (i * m)] = 0.0;

                }
            constraintRHS[i] = 0.0;

         } //fill each array with zeros

   double sigma = data[0];
   double Xo = data[1];
   double Yo = data[2];

    //Now want constraints sigma = data[0], Xo = data[1] Yo = data[3]. constraintMatrix is 3x5
    //So only interested in elements (1,2), (2,4) and (3,5) - 2nd, 9th and 15th elements (levmar uses row major order)
    constraintMatrix[1] = 1;
    constraintRHS[0] = sigma; //sigma = input

    constraintMatrix[8] = 1;
    constraintRHS[1] = Xo; //Xo = input

    constraintMatrix[14] = 1;
    constraintRHS[2] = Yo; //Yo = input


    }

///==================ELLIPTICAL CONSTRAINT FUNCTIONS==================


//basic bounding constraints for an elliptical gaussian, bounds for sigma and position specified by the user
void constraints_basic_elliptical(double *lb, double *ub, int numBounds, int numEqualityConstraints, double *weights, int m, int framesize, double *data)
{


    /// bounds - A,b,theta hardcoded, others specified from Matlab input

    //A constraints  - A is greater than or equal to zero (non-negative)
    lb[0] = 0.0;
    ub[0] = 65535.0;


    //sigma constraints
    //sigma_x
    lb[1] = data[numEqualityConstraints]; //equality constraints are loaded into data first, so the bounds come after
    ub[1] = data[numEqualityConstraints+1];

    //sigma_y
    lb[2] = data[numEqualityConstraints+2]; //equality constraints are loaded into data first, so the bounds come after
    ub[2] = data[numEqualityConstraints+3];


    //b constraints - b is any value
    lb[3] = -65535.0;
    ub[3] = 65535.0;

    //Xo and Yo constraints
    lb[4] = data[numEqualityConstraints+4];
    ub[4] = data[numEqualityConstraints+5]; //Xo

    lb[5] = data[numEqualityConstraints+6];
    ub[5] = data[numEqualityConstraints+7]; //Yo

    //theta constraints - 0 to 2pi


    const double pi = 3.14159265358979323846;

    //lb[6] = 0.0;
    //ub[6] = 2*pi;
    lb[6] = -2*pi;
    ub[6] = 2*pi; // add a bit more wiggle room
    //lb[6] = -50*pi;
    //ub[6] = 50*pi; // add a bit more wiggle room

    return;
    }

///FIXED Xo, Yo (position) (elliptical)
void constraints_fixedposition_elliptical(double *constraintMatrix, double *constraintRHS, int numBounds, int numEqualityConstraints, int m, double* data)

    {

        for (int i = 0; i < numEqualityConstraints; i++)  //fill constraintMatrix and RHS row by row
         {

             for (int j = 0; j < m; j++)
                {

                 constraintMatrix[j + (i * m)] = 0.0;

                }
            constraintRHS[i] = 0.0;

         } //fill each array with zeros

   double Xo = data[0];
   double Yo = data[1];
                               //(first 0->numBounds-1 of data contain the bounding constraints
    //Now want constraints Xo = data[numBounds] Yo = data[numBounds+1]. constraintMatrix is 2x7
    //So only interested in elements (1,5) and (2,6) - 5th and 13th elements (levmar uses row major order)
    constraintMatrix[4] = 1;
    constraintRHS[0] = Xo; //Xo = input

    constraintMatrix[12] = 1;
    constraintRHS[1] = Yo; //Yo = input

    return;
    }

//write new functions here


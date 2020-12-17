//fit functions for contiunuous circular gaussian, continuous elliptical gaussian and background only models
//each includes a function to calculate their analytic jacobian
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



///==============CIRCULAR GAUSSIAN FITTING=============================

     void
     gauss_f (double *p, double *x, int m, int n, void *data )
     {
       double ndbl = n;
       double framesize = sqrt(ndbl);
       //gets data
       int xsize = framesize; //n should be a square number so this is ok
       int ysize = xsize;  //assumes square data area


    //grabs the current guesses
       int p_count = 0;
       double A = p[p_count++];
       double sigma_xy = p[p_count++];
       double b = p[p_count++];
       double Xo = p[p_count++];
       double Yo = p[p_count++];

    //initialises other useful variables
       double X, Y, e;
       int element;

    //loops over all pixels in the image
       for(int y_coord = 0; y_coord < ysize; y_coord++)
         {
       for(int x_coord = 0; x_coord < xsize; x_coord++)
         {
           X = x_coord + 1; //+1 because we start at 0 in the for loop
                            //but the real data will be from a matrix indexed to start at [1,1]
           Y = y_coord + 1;

           element = x_coord + (y_coord * ysize);
           //y_coord *y_size so it starts at 0, then 1*ysize...so element[2,1] in the original matrix is
           //element y_size+1 in f etc

           //=======Model Function - amplitude * 2D gaussian + background==========
             e = exp (-((pow(X-Xo,2) + pow(Y-Yo,2)) * 0.5 * pow(sigma_xy,-2))); //2D gaussian
             x[element] = (A * e) + b;
           //======================================================================
         }
         }



     }

     void
     gauss_jac (double *p, double *jac, int m, int n, void *data)
      //calculates the Jacobian vector of derivatives for each of the variables to be fitted

     {
       double ndbl = n;
       double framesize = sqrt(ndbl);
       //gets data
       int xsize = framesize; //n should be a square number so this is ok
       int ysize = xsize;  //assumes square data area




    //grabs the initial guesses
       int p_count = 0;
       double A = p[p_count++];
       double sigma_xy = p[p_count++];
       double b = p[p_count++];
       double Xo = p[p_count++];
       double Yo = p[p_count++];
       //initialises other useful variables
       double X, Y, e;
       int element;
       int j = 0;

    //loops over all pixels in the image
       for (int y_coord = 0; y_coord < ysize; y_coord++)
         {
       for (int x_coord = 0; x_coord < xsize; x_coord++)
         {
           // Jacobian matrix J(i,j) = dfi / dxj,
           // fi = 2d gaussian with background
           // and the xj are the parameters (A,sigma,b,Xo,Yo)

           X = x_coord + 1;
           Y = y_coord + 1;
           element = x_coord + (y_coord * ysize);
           e = exp (-((pow(X-Xo,2) + pow(Y-Yo,2)) * 0.5 * pow(sigma_xy,-2)));



           jac[j++] = e; //dF/dA
           ///can put if variant =! 2 then run this next line style thing
           jac[j++] = pow(sigma_xy,-3) * ( pow(X-Xo,2) + pow(Y-Yo,2) ) * A * e;  //dF/d(sigma)
           jac[j++] = 1.0;                                    //dF/db
           jac[j++] = pow(sigma_xy,-2) * (X-Xo) * A * e; //dF/dXo
           jac[j++] = pow(sigma_xy,-2) * (Y-Yo) * A * e; //dF/dYo
         }
         }

     }

///------------------------------------------------------///
///                 ELLIPTICAL FITTING                   ///

     void
     elliptical_gauss_f (double *p, double *x, int m, int n, void *data )
     //calculates the values of the model equation to be fitted and the residuals on each point
     //for an elliptical gaussian (free params A, sigma_x, sigma_y, b, Xo, Yo, theta

     {

       double ndbl = n;
       double framesize = sqrt(ndbl);
       //gets data
       int xsize = framesize; //n should be a square number so this is ok
       int ysize = xsize;  //assumes square data area


    //grabs the current guesses
       int p_count = 0;
       double A = p[p_count++];
       double sigma_x = p[p_count++];
       double sigma_y = p[p_count++];
       double b = p[p_count++];
       double Xo = p[p_count++];
       double Yo = p[p_count++];
       double theta = p[p_count++];

      //initialises variables to calculate model function
       double X, xprime, Y, yprime, e, xdenom, ydenom; //useful subparts of the elliptical Gaussian
       int element;

       xdenom = sqrt(2.0) * sigma_x;
       ydenom = sqrt(2.0) * sigma_y;

    //loops over all pixels in the image
       for(int y_coord = 0; y_coord < ysize; y_coord++)
         {
       for(int x_coord = 0; x_coord < xsize; x_coord++)
         {
           // Model Yi = 2D gaussian with background and std dev sigma_xy
           X = x_coord + 1; //+1 because we start at 0 in the for loop but the real data will be from a matrix indexed to start at [1,1]
           Y = y_coord + 1;

           element = x_coord + (y_coord * ysize);
           //y_coord *y_size so it starts at 0, then 1*ysize...so element[2,1] in the original matrix is
           //element y_size+1 in f etc

        //calculates useful sub-parts of the elliptical gaussian
           xprime = (X-Xo)*cos(theta) - (Y-Yo)*sin(theta);
           yprime = (X-Xo)*sin(theta) + (Y-Yo)*cos(theta);

           e = exp((-(pow(xprime/xdenom,2)))-(pow(yprime/ydenom,2)));


           x[element] = (A * e) + b;
         }
         }


       return;
     }


     void
     elliptical_gauss_jac (double *p, double *jac, int m, int n, void *data)
      //calculates the Jacobian vector of derivatives for each of the variables to be fitted

     {
       double ndbl = n;
       double framesize = sqrt(ndbl);
       //gets data
       int xsize = framesize; //n should be a square number so this is ok
       int ysize = xsize;  //assumes square data area


    //grabs the initial guesses
       int p_count = 0;
       double A = p[p_count++];
       double sigma_x = p[p_count++];
       double sigma_y = p[p_count++];
       double b = p[p_count++];
       double Xo = p[p_count++];
       double Yo = p[p_count++];
       double theta = p[p_count++];
//	cout << A << " " << sigma_x << " " << sigma_y << endl;

       //initialises other useful variables
       int j = 0;
       double X, xprime, Y, yprime, e, xdenom, ydenom; //useful subparts of the elliptical Gaussian
       int element;

       xdenom = sqrt(2.0) * sigma_x;
       ydenom = sqrt(2.0) * sigma_y;

    //loops over all pixels in the image
       for (int y_coord = 0; y_coord < ysize; y_coord++)
         {
       for (int x_coord = 0; x_coord < xsize; x_coord++)
         {
           // Jacobian matrix J(i,j) = dfi / dxj,
           // where fi = (Yi - yi)/
           // Yi = 2d gaussian with background
           // and the xj are the parameters (A,sigma,b,Xo,Yo)

           X = x_coord + 1;
           Y = y_coord + 1;
           element = x_coord + (y_coord * ysize);


        //calculates useful sub-parts of the elliptical gaussian
           xprime = (X-Xo)*cos(theta) - (Y-Yo)*sin(theta);
           yprime = (X-Xo)*sin(theta) + (Y-Yo)*cos(theta);

           e = exp((-(pow(xprime/xdenom,2)))-(pow(yprime/ydenom,2)));


           jac[j++] = e; //dF/dA

           jac[j++] =  (A * e) * (pow(xprime,2) * pow(sigma_x,-3)); //dF/dsigma_x
           jac[j++] =  (A * e) * (pow(yprime,2) * pow(sigma_y,-3)); //dF/dsigma_y

           jac[j++] = 1; //dF/db

           jac[j++] = (A * e) * ( (xprime*cos(theta)*pow(sigma_x,-2)) + (yprime*sin(theta)*pow(sigma_y,-2)) ); //dF/dXo
           jac[j++] = (A * e) * ( (yprime*cos(theta)*pow(sigma_y,-2)) - (xprime*sin(theta)*pow(sigma_x,-2)) ); //dF/dYo

           jac[j++] = (-A * e) *( (-xprime * pow(sigma_x,-2)) * ((X-Xo)*sin(theta) + (Y-Yo)*cos(theta)) + (yprime*pow(sigma_y,-2))*((X-Xo)*cos(theta) - (Y-Yo)*sin(theta)) ); //dF/dtheta
         }
         }
       return;
     }


///------------------------------------------------------///
///                 BACKGROUND FITTING                   ///

        void
        back_f (double *p, double *x, int m, int n, void *data )      //calculates the residuals on each point using the model F = constant

     {
       double ndbl = n;
       double framesize = sqrt(ndbl);
       //gets data
       int xsize = framesize; //n should be a square number so this is ok
       int ysize = xsize;  //assumes square data area
      // double *data_array = ((struct data *)data)->data_array;


       //only need b for the model
       double b = p[0]; //zero here because we make a separate array to hold b alone
       int element;

       for(int y_coord = 0; y_coord < ysize; y_coord++)
         {
       for(int x_coord = 0; x_coord < xsize; x_coord++)

         {
           // Model x[i] =  const
           element = x_coord + (y_coord * ysize);
           x[element] =  b;
         }
         }

     }

     void
     back_jac (double *p, double *jac, int m, int n, void *data)       //calculates the Jacobian matrix of derivatives for each of the variables to be fitted
     {

       double ndbl = n;
       double framesize = sqrt(ndbl);
       //gets data
       int xsize = framesize; //n should be a square number so using an int here is ok
       int ysize = xsize;  //assumes square data area
      // double *data_array = ((struct data *)data)->data_array;



       int j = 0;

       for (int y_coord = 0; y_coord < ysize; y_coord++)
         {
       for (int x_coord = 0; x_coord < xsize; x_coord++)
         {
            //Jacobian vector dx[i]/dp[j] over all i, j.

           jac[j++] = 1; //jacobian is not very interesting

         }
         }
     }


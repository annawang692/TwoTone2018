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
%GNU General Public License (below) for more details.
%
*/
//Code for generating an initial guess of parameters for a continuous gaussian model using data provided

//Uses mostly same naming convention as main program but array for storing measured data, x[n],
//is now called data_array[n] for clarity due to existence of x,y coords


    void initialfit(char *variant, double p[], double data_array[], int n, double *data )
    //writes the initial guess fitting variables to array p[m] (and calculates them from data)
    //n is no of observations, framesize the width of the data area (n = framesize^2) - assume square box


    //Currently have program set up so that based on 'variant' it knows which set of variables you are using:
    //free ,fixedwidth, fixedposition, fixedwidthposition

    //for fixed conditions, data contains all the parameters needed following general convention of having
    //any array containing some or all of the independent variables in this order - amplitude, sigma, background, position
    //with x being in front of y
    //so for fixed width position data[0] = sigma data[1] = Xo data[2] = Yo
    //p[m] has same conventions

    //If you wanted to fix/unfix something then hack/write a new function, add it at the top
    //and include another condition and function definition below

    //all find_ functions assume circular gaussian model, all ellipt_find_ assume elliptical gaussian model


        {
    //self explanatorily named  guessing functions - can see what they depend on already having calculated from arguments
        void find_AbXY (double p[], double data_array[], int framesize );

        void find_Ab_positionknown(double p[], double data_array[], int framesize, double Xo, double Yo);

        void find_b (double p[], double data_array, int framesize);

        void find_sigma (double p[], double data_array[], int framesize, double Ao, double b, double Xo, double Yo);

        void find_AbXYtheta_elliptical (double p[], double data_array[], int framesize);

        void find_Abtheta_positionknown_elliptical(double p[], double data_array[], int framesize, double Xo, double Yo);

        void find_sigma_elliptical (double p[], double data_array[], int framesize, double Ao, double b, double Xo, double Yo);


        double ndbl = n;               //bit hacky - matlab won't square root an int so convert int to double
        int framesize = floor(sqrt(ndbl)); //then sqrt it and then convert back to int
        double Ao, b, Xo, Yo; // sigma_xy, sigma_x, sigma_y, theta;


        if(variant == "free")
            {
                find_AbXY (p, data_array, framesize);
                find_sigma(p, data_array, framesize, p[0], p[1], p[3], p[4]);
            }

        else if(variant == "fixedposition" )
            {
                Xo = data[0];
                Yo = data[1];
                p[3] = Xo; //Xo
                p[4] = Yo; //Yo

                find_Ab_positionknown(p, data_array, framesize, Xo, Yo);
                Ao = p[0];
                b = p[1];

                find_sigma(p, data_array, framesize, Ao, b, Xo, Yo);
            }

        else if(variant == "fixedsigma" )
            {
                p[1] = data[0]; //sigma
                find_AbXY (p, data_array, framesize);

            }

        else if(variant == "fixedsigmaposition" )
            {
                p[1] = data[0]; //sigma
                p[3] = data[1]; //Xo
                p[4] = data[2]; //Yo

                find_Ab_positionknown(p, data_array, framesize, p[3], p[4]);
            }

        else if(variant == "free_elliptical")
            {
                find_AbXYtheta_elliptical (p, data_array, framesize);
                find_sigma_elliptical (p, data_array, framesize, p[0], p[3], p[4], p[5]);
            }

        else if(variant == "fixedposition_elliptical")
            {
                Xo = data[0];
                Yo = data[1];
                p[4] = Xo;
                p[5] = Yo;

                find_Abtheta_positionknown_elliptical (p, data_array, framesize, Xo, Yo);
                find_sigma_elliptical (p, data_array, framesize, p[0], p[3], p[4], p[5]);
            }

        else
            {mexErrMsgTxt("initial fit failure - given variant not found");}


        }



//subfunctions here are called to determine initial fit variables

    void find_AbXY (double p[], double data_array[], int framesize )
   //Returns Ao, b, Xo and Yo in a 4 element array, in that order
    {

        double Ao;
        double A_raw = data_array[0];
        double b_sum = 0.0;
        double Xo = 0;
        double Yo = 0;


        //first we want to find the background, b (p[3]) and the raw max amplitude
        //we do this by guessing that the background will be close in value to the mean value of intensity
        //and the raw max amplitude will be the max valued datapoint

    for(int y = 0; y < framesize; y++)
        {
    for(int x = 0; x < framesize; x++) //start at 1 because we already set A_raw, b, to y[0]
        {
            if( data_array[((y*framesize) + x)]> A_raw ) //evaluates whether the new datapoint has a bigger value than the previous largest, if so, sets X_POS to its index and A_raw to its value
            {
                Xo = x+1;
                Yo = y+1;
                A_raw = data_array[((y*framesize) + x)];
            }

        b_sum = b_sum + data_array[((y*framesize) + x)];  //evaluates whether the new datapoint is a new minimum and if so sets b to it

        }
        }
    //Now we guess Ao the adjusted amplitude will be ~A_raw-b
    double b = b_sum/(framesize*framesize); //get avg data point value
    Ao = A_raw - b;

    //now we have Ao, sigma/, b and X_POS and Y_POS
    //so return them as array

    p[0] = Ao;
    p[2] = b;
    p[3] = Xo;
    p[4] = Yo;

    return;


}


    void find_Ab_positionknown(double p[], double data_array[], int framesize, double Xo, double Yo)

    {
        double A_raw = data_array[0];
        double b_sum = 0.0;

        for(int y = 0; y < framesize; y++)
         {
        for(int x = 0; x < framesize; x++) //start at 1 because we already set A_raw, b, to y[0]
         {

             if( data_array[((y*framesize) + x)]> A_raw ) //evaluates whether the new datapoint has a bigger value than the previous largest, if so, sets X_POS to its index and A_raw to its value
              {
              A_raw = data_array[((y*framesize) + x)];
              }

         b_sum = b_sum + data_array[((y*framesize) + x)];  //evaluates whether the new datapoint is a new minimum and if so sets b to it

         }
         }

    //Now we guess Ao the adjusted amplitude will be ~A_raw-b
        double b = b_sum/(framesize*framesize); //get avg data point value - divide by num of pixels
        double Ao = A_raw - b;

        p[0] = Ao;
        p[2] = b;

        return;



    }


    void find_b (double p[], double data_array[], int framesize)

    {

    double b_sum = 0.0;

        for(int y = 0; y < framesize; y++)
         {
        for(int x = 0; x < framesize; x++) //start at 1 because we already set A_raw, b, to y[0]
         {

         b_sum = b_sum + data_array[((y*framesize) + x)];  //evaluates whether the new datapoint is a new minimum and if so sets b to it

         }
         }

    //Now we guess Ao the adjusted amplitude will be ~A_raw-b
    double b = b_sum/(framesize*framesize); //get avg data point value

    p[2] = b;

    return;

    }



    void find_sigma (double p[], double data_array[], int framesize, double Ao, double b, double Xo, double Yo)
    //guesses at sigma using other 4 parameters - returns sigma

    {

        double minus_one = -1;
        double expdecay = exp(minus_one);

        int searchside = 1; //sets up the search variable to search the larger area (in case it is cut off on one side)
        if(Xo > framesize/2)
        { searchside = -1; }


        double Fx;
        int j = Xo;

        do {

            j = j + searchside;   //increments the search
            Fx = data_array[j];           //sets Fx to the datavalue at the current datapoint

                    if(j > framesize || j<0)  //error check in case
                        {
                          j = Xo+1;
                          break;
                        }

          } while(Fx > ((Ao * expdecay) + b)); //sigma should be distance from peak at which F = b + A * exp(-1)

        double sigma_x = sqrt( sqrt (pow(((Xo - j)*0.5),2))); //finds X_POS - j, which is 2 sigma squared,halves it, takes the root of the modulus of it to get sigma

        searchside = 1; //sets up the search variable to search the larger area (in case it is cut off on one side)
        if(Yo > framesize/2)
        { searchside = -1; }


        double Fy;
        j = Yo;

        do {

            j = j + searchside;   //increments the search
            Fy = data_array[j];           //sets Fy to the datavalue at the current datapoint

                    if(j > framesize || j<0)  //error check in case j has gone outside data area
                        {
                          j = Yo-1;
                          break;
                        }

          } while(Fy > ((Ao * expdecay) + b));

        double sigma_y = sqrt( sqrt (pow(((Yo - j)*0.5),2))); //finds X_POS - j, which is 2 sigma squared,halves it, takes the root of the modulus of it to get sigma

        double sigma_xy = (sigma_x + sigma_y)/2; //takes avg

        p[1] = sigma_xy;

        return;

    }


    void find_AbXYtheta_elliptical (double p[], double data_array[], int framesize)
        //doesn't actually do anything different to circular one as just sets theta to 0
        //but updates the p[] indices at the end to account for having 2 values of sigma for x and y
    {


        double Ao,theta, b;
        int Xo, Yo;

        //Let theta = 0 - assume an x-y centered eliipse to begin with
        theta = 0;

        double A_raw = data_array[0];
        double b_sum = 0.0;

        b = data_array[0];

        Xo = 0;
        Yo = 0;

        //first we want to find the background, b (x_init[3]) and the raw max amplitude
        //we do this by guessing that the background will be close in value to the mean value of intensity
        //and the raw max amplitude will be the max valued datapoint

    for(int y = 0; y < framesize; y++)
    {
    for(int x = 0; x < framesize; x++) //start at 1 because we already set A_raw, b, to y[0]
    {
      if( data_array[((y*framesize) + x)]> A_raw ) //evaluates whether the new datapoint has a bigger value than the previous largest, if so, sets X_POS to its index and A_raw to its value
       {
          Xo = x+1;
          Yo = y+1;
          A_raw = data_array[((y*framesize) + x)];
       }

        b_sum = b_sum + data_array[((y*framesize) + x)];  //evaluates whether the new datapoint is a new minimum and if so sets b to it

        }
        }
        //Now we guess Ao the adjusted amplitude will be ~A_raw-b
        b = b_sum/(framesize*framesize); //get avg data point value
        Ao = A_raw - b;

        //to find sigma we want to search on the side with more space (quickest to search only 1 side although searching both sides should be easy)

        //now we have Ao, b, Xo, Yo, theta
        //so put them back into the array

        p[0] = Ao;
        p[3] = b;
        p[4] = Xo;
        p[5] = Yo;
        p[6] = theta;

        return;
}
    void find_Abtheta_positionknown_elliptical(double p[], double data_array[], int framesize, double Xo, double Yo)
    //copy of circular version with indices updated on output to account for sigma_x and sigma_y and theta set to 0
        {
        double A_raw = data_array[0];
        double b_sum = 0.0;
        double theta = 0.0;

        for(int y = 0; y < framesize; y++)
         {
        for(int x = 0; x < framesize; x++) //start at 1 because we already set A_raw, b, to y[0]
         {

             if( data_array[((y*framesize) + x)]> A_raw ) //evaluates whether the new datapoint has a bigger value than the previous largest, if so, sets X_POS to its index and A_raw to its value
              {
              A_raw = data_array[((y*framesize) + x)];
              }

         b_sum = b_sum + data_array[((y*framesize) + x)];  //evaluates whether the new datapoint is a new minimum and if so sets b to it

         }
         }

    //Now we guess Ao the adjusted amplitude will be ~A_raw-b
        double b = b_sum/(framesize*framesize); //get avg data point value - divide by num of pixels
        double Ao = A_raw - b;

        p[0] = Ao;
        p[3] = b;
        p[6] = theta;

        return;

    }

    void find_sigma_elliptical (double p[], double data_array[], int framesize, double Ao, double b, double Xo, double Yo)

    {

        double minus_one = -1;
        double expdecay = exp(minus_one);

        int searchside = 1; //sets up the search variable to search the larger area (in case it is cut off on one side)
        if(Xo > framesize/2)
        { searchside = -1; }


        double Fx;
        int j = Xo;

        do {

            j = j + searchside;   //increments the search
            Fx = data_array[j];           //sets Fx to the datavalue at the current datapoint

                    if(j > framesize || j<0)  //error check in case
                        {
                          j = Xo+1;
                          break;
                        }

          } while(Fx > ((Ao * expdecay) + b)); //sigma should be distance from peak at which F = b + A * exp(-1)

        double sigma_x = sqrt( sqrt (pow(((Xo - j)*0.5),2))); //finds Xo - j, which is 2 sigma squared,halves it, takes the root of the modulus of it to get sigma

        searchside = 1; //sets up the search variable to search the larger area (in case it is cut off on one side)
        if(Yo > framesize/2)
        { searchside = -1; }


        double Fy;
        j = Yo;

        do {

            j = j + searchside;   //increments the search
            Fy = data_array[j];           //sets Fy to the datavalue at the current datapoint

                    if(j > framesize || j<0)  //error check in case j has gone outside data area
                        {
                          j = Yo-1;
                          break;
                        }

          } while(Fy > ((Ao * expdecay) + b));

        double sigma_y = sqrt( sqrt (pow(((Yo - j)*0.5),2))); //finds Yo - j, which is 2 sigma squared,halves it, takes the root of the modulus of it to get sigma

        p[1] = sigma_x;
        p[2] = sigma_y;

        return;
    }






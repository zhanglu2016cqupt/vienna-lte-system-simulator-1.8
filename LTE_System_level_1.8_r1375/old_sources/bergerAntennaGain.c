/* file: bergerAntennaGain.c

   The calling syntax is:

      [output] = bergerAntennaGain(theta, mean_antenna_gain)

      output = antenna gain

   Josep Colom Ikuno, INTHFT
   josep.colom@nt.tuwien.ac.at
   www.nt.tuwien.ac.at

*/

#include <mex.h>
#include <math.h>
#include <stdlib.h>

/* Input Arguments */
#define INPUT1       prhs[0]
#define INPUT2       prhs[1]

/* Output Arguments */
#define OUTPUT       plhs[0]

/* main function that interfaces with MATLAB */
void mexFunction(
				 int            nlhs,
				 mxArray       *plhs[],
				 int            nrhs,
				 const mxArray *prhs[] )
{
    double    theta;
    double    mean_antenna_gain;
    double    *gain;
    double     temp;
    
    /* Check if the input is of the correct type */
    if ( (mxIsDouble(INPUT1) != 1) && (mxIsDouble(INPUT2) != 1)) {
        mexErrMsgTxt("Input must be double.");
    }
    else {
        /* first input is the data word */
        theta = *mxGetPr(INPUT1);
        mean_antenna_gain = *mxGetPr(INPUT2);
        
        /* create the output vector */
        OUTPUT = mxCreateDoubleMatrix(1,1,mxREAL);
        gain = mxGetPr(OUTPUT);
        
        /* Populate the output */
        temp = 12*pow((theta/70),2);
        temp = temp<20?temp:20;
        *gain = -temp + mean_antenna_gain;
    }
    
    return;
}

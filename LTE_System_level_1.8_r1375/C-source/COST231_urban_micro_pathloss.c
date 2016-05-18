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
#define INPUT3       prhs[2]
#define INPUT4       prhs[3]
#define INPUT5       prhs[4]
#define INPUT6       prhs[5]
#define INPUT7       prhs[6]
#define INPUT8       prhs[7]

/* Output Arguments */
#define OUTPUT       plhs[0]

/* main function that interfaces with MATLAB */
void mexFunction(
				 int            nlhs,
				 mxArray       *plhs[],
				 int            nrhs,
				 const mxArray *prhs[] )
{
    double    *distance;
    double    actual_distance;
    double    frequency;
    double    L_0;
    double    h_base, h_roof, h_mobile;
    double    delta_hbase, delta_hmobile;
    double    phi, w, b;
    double    L_ori, L_rts, L_bsh, L_msd;
    double    k_a, k_d, k_f;
    double    *pl_NLOS;
    mwSize    input_length;
    mwSize    i;
    mwSize    *input_dimensions;
    mwSize    number_of_input_dimensions;
    
    double    temp;
    
    /* No checking here, BE CAREFUL! */
    distance  =  mxGetPr(INPUT1);
    frequency = *mxGetPr(INPUT2);
    h_roof    = *mxGetPr(INPUT3);
    h_base    = *mxGetPr(INPUT4);
    h_mobile  = *mxGetPr(INPUT5);
    phi       = *mxGetPr(INPUT6);
    w         = *mxGetPr(INPUT7);
    b         = *mxGetPr(INPUT8);
    input_length = mxGetNumberOfElements(INPUT1);
    input_dimensions = mxGetDimensions(INPUT1);
    number_of_input_dimensions = mxGetNumberOfDimensions(INPUT1);
    
    /* create the output vector */
    OUTPUT  = mxCreateDoubleMatrix(1,input_length,mxREAL);
    OUTPUT =  mxCreateNumericArray(number_of_input_dimensions,input_dimensions,mxDOUBLE_CLASS,mxREAL);
    pl_NLOS = mxGetPr(OUTPUT);
    
    // Calculations are done in freq in MHz
    frequency = frequency/1000000;
    
    /* Populate the output */
    for(i=0;i<input_length;i++) {
        // Calculations are done in Km
        actual_distance = (double)(distance[i]/1000);
        
        // NLOS component
        L_0 = 32.4 + 20*log10(actual_distance) + 20*log10(frequency); // free space loss
        
        delta_hmobile = h_roof - h_mobile;
        if ((0 <= phi) && (phi < 35)) {
            L_ori = -10 + 0.354*phi;
        } else if((35 <= phi) && (phi < 55)) {
            L_ori = 2.5 + 0.075*(phi - 35);
        } else if ((55 <= phi) && (phi < 90)) {
            L_ori = 4.0 - 0.114*(phi - 55);
        } else if ((55 <= phi) && (phi < 90)) {
            L_ori = 4.0 - 0.114*(phi - 55);
        } else {
            mexErrMsgTxt("Wrong wave propagation angle.");
        }
        
        L_rts = -16.9 - 10*log10(w) + 10*log10(frequency) + 20*log10(delta_hmobile) + L_ori; // roof-top-to-street diffracton and scatter loss
        
        delta_hbase = h_base - h_roof;
        if (h_base > h_roof) {
            L_bsh = -18*log10(1 + delta_hbase);
        } else {
            L_bsh = 0;
        }
        
        if(h_base > h_roof) {
            k_a = 54;
        } else if((actual_distance >= 0.5) && (h_base <= h_roof)) {
            k_a = 54 - 0.8*delta_hbase;
        } else if((actual_distance < 0.5) && (h_base <= h_roof)) {
            k_a = 54 - 0.8*delta_hbase*actual_distance/0.5;
        }
        
        if(h_base > h_roof) {
            k_d = 18;
        } else {
            k_d = 18 - 15*delta_hbase/h_roof;
        }
        
        k_f = -4 + 1.5*(frequency/925 - 1); // for metropolitan areas
        // for the other two, the COST231-WI model does not apply ;)
        L_msd = L_bsh + k_a + k_d*log10(actual_distance) + k_f*log10(frequency) - 9*log10(b); //multiple screen diffraction loss
        
        temp = (L_rts + L_msd)>0 ? L_rts + L_msd:0;
        pl_NLOS[i] = L_0 + temp; // NLOS pathloss
    }
    
    return;
}

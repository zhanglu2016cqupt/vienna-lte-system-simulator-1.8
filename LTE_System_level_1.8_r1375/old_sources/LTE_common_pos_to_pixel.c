/* file: LTE_common_pos_to_pixel.c. Translates from an absolut position to a pixel position

   The calling syntax is:

      [ pos_pixel pos_pixel_exact] = LTE_common_pos_to_pixel( pos, roi_min, data_res)

% Converts a position in absolute values to a pixel position
% (c) Josep Colom Ikuno, INTHFT, 2008
% input:    pos             ... [x,y] Position to convert
%           roi_min         ... [x,y] Lower-leftmost corner of the ROI
%           data_res        ... meters/pixel, resolution of the map
% output:   pos_pixel       ... [x,y] Pixel position (1-indexed)
%           pos_pixel_exact ... [x,y] Pixel position, exact value. Can be
%                               used to perform interpolation. As with
%                               pos_pixel, it is one-indexed.

*/

#include <mex.h>
#include <math.h>
#include <stdlib.h>

/* Input Arguments */
#define INPUT1       prhs[0]
#define INPUT2       prhs[1]
#define INPUT3       prhs[2]

/* Output Arguments (depending on nlhs) */
#define OUTPUT1       plhs[0]
#define OUTPUT2       plhs[1]

/* main function that interfaces with MATLAB */
void mexFunction(
				 int            nlhs,
				 mxArray       *plhs[],
				 int            nrhs,
				 const mxArray *prhs[] )
{
    double    *pos;
    double    *roi_min;
    double    data_res;
    double    *pos_pixel;
	double    *pos_pixel_exact;
    
    /* No type checking. Remember that input MUST be double!! */
    
    pos      = mxGetPr(INPUT1);
	roi_min  = mxGetPr(INPUT2);
	data_res = *mxGetPr(INPUT3);
        
        /* create the output vector */
        OUTPUT1 = mxCreateDoubleMatrix(1,2,mxREAL);
        pos_pixel = mxGetPr(OUTPUT1);
		
		pos_pixel[0] = floor((pos[0]-roi_min[0])/data_res)+1;	
		pos_pixel[1] = floor((pos[1]-roi_min[1])/data_res)+1;
		
		if(nlhs>1) {
			OUTPUT2 = mxCreateDoubleMatrix(1,2,mxREAL);
			pos_pixel_exact = mxGetPr(OUTPUT2);
			
			pos_pixel_exact[0] = ((pos[0]-roi_min[0])/data_res)+1;	
			pos_pixel_exact[1] = ((pos[1]-roi_min[1])/data_res)+1;
		}
    
    return;
}

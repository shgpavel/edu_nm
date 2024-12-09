#include <math.h>

__global__ void rk2wc_cuda(double t0, double t_end, double h,
													 double tol, double A, double B,
													 double a21, double b1, double b2,
													 double *y_initial, double *y_values) {
		
	int system_idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	double *y_current = y_initial + 2 * system_idx;
	double *y_all_steps = y_values + system_idx;
	
	double y1 = y_current[0];
	double y2 = y_current[1];

	y_all_steps[0] = y1;
	y_all_steps[num_steps * 2 * gridDim.x] = y2;

	double h_2 = h * 0.5f;
	
	for (unsigned step = 1; step < num_steps; ++step) {
		/* f */
		double f_y1 = A * y2;
		double f_y2 = -B * y1;
		
		/* associative part */
		double ap_1 = a21 * h * f_y1;
		double ap_2 = a21 * h * f_y2;
		
		/* first part of full step */
		double y_1p1 = y1 + ap_1;
		double y_2p1 = y2 + ap_2;
		
		/* apply f 2nd time */
		double f2_y1 = A * y_1p1;
		double f2_y2 = -B * y_2p1;
    
		/* first part of full step */
		double y_1p1_2 = y1 + ap_1 / 2;
		double y_2p1_2 = y2 + ap_2 / 2;
		
		/* apply f 2nd time */
		double f2_y1_2 = A * y_1p1_2;
		double f2_y2_2 = -B * y_2p1_2;
		
		double y1_full = y1 + b1 * h * f_y1 + b2 * h * f2_y1;
		double y2_full = y2 + b1 * h * f_y2 + b2 * h * f2_y2;
		
		double y1_half = y1 + b1 * h_2 * f_y1 + b2 * h_2 * f2_y1_2;
		double y2_half = y2 + b1 * h_2 * f_y1 + b2 * h_2 * f2_y2_2;
		
		double error1 = fabsf(y1_full - y1_half);
		double error2 = fabsf(y2_full - y2_half);
		double error = fmaxf(error1, error2) / 3.0f;
		
		y_all_steps[step * 2 * gridDim.x] = y1;
		y_all_steps[step * 2 * gridDim.x + 1] = y2; 

		/*
		if (error > tol) {
			h *= 0.5f;
		} else {
			
		}
		*/
	}
}

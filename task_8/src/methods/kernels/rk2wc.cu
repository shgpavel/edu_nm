#include <math.h>

__global__ void rk2wc_cuda(float h, float A, float B,
													 float a21, float b1, float b2,
													 float tol, unsigned num_steps,
													 float *y_initial, float *y_values) {
	
	int system_idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	float *y_current = y_initial + 2 * system_idx;
	float *y_all_steps = y_values + system_idx;
	
	float y1 = y_current[0];
	float y2 = y_current[1];

	y_all_steps[0] = y1;
	y_all_steps[num_steps * 2 * gridDim.x] = y2;

	float h_2 = h * 0.5f;
	
	for (unsigned step = 1; step < num_steps; ++step) {
		/* f */
		float f_y1 = A * y2;
		float f_y2 = -B * y1;
		
		/* associative part */
		float ap_1 = a21 * h * f_y1;
		float ap_2 = a21 * h * f_y2;
		
		/* first part of full step */
		float y_1p1 = y1 + ap_1;
		float y_2p1 = y2 + ap_2;
		
		/* apply f 2nd time */
		float f2_y1 = A * y_1p1;
		float f2_y2 = -B * y_2p1;
    
		/* first part of full step */
		float y_1p1_2 = y1 + ap_1 / 2;
		float y_2p1_2 = y2 + ap_2 / 2;
		
		/* apply f 2nd time */
		float f2_y1_2 = A * y_1p1_2;
		float f2_y2_2 = -B * y_2p1_2;
		
		float y1_full = y1 + b1 * h * f_y1 + b2 * h * f2_y1;
		float y2_full = y2 + b1 * h * f_y2 + b2 * h * f2_y2;
		
		float y1_half = y1 + b1 * h_2 * f_y1 + b2 * h_2 * f2_y1_2;
		float y2_half = y2 + b1 * h_2 * f_y1 + b2 * h_2 * f2_y2_2;
		
		float error1 = fabsf(y1_full - y1_half);
		float error2 = fabsf(y2_full - y2_half);
		float error = fmaxf(error1, error2) / 3.0f;
		
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

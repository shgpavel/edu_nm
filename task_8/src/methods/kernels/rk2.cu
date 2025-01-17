#include <math.h>

__global__ void rk2_cuda(float h, float A, float B,
												 float a21, float b1, float b2,
												 unsigned num_steps, float *y_initial,
												 float *y_values) {
    
	int system_idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	float *y_current = y_initial + 2 * system_idx;
	float *y_all_steps = y_values + system_idx;
	
	float y1 = y_current[0];
	float y2 = y_current[1];

	y_all_steps[0] = y1;
	y_all_steps[num_steps * 2 * gridDim.x] = y2;
	
	for (unsigned step = 1; step < num_steps; ++step) {
		float k1_y1 = A * y2;
		float k1_y2 = -B * y1;
    
		float y_temp1 = y1 + a21 * h * k1_y1;
		float y_temp2 = y2 + a21 * h * k1_y2;
		
		float k2_y1 = A * y_temp2;
		float k2_y2 = -B * y_temp1;
    
		y1 = y1 + h * (b1 * k1_y1 + b2 * k2_y1);
		y2 = y2 + h * (b1 * k1_y2 + b2 * k2_y2);
    
		y_all_steps[step * 2 * gridDim.x] = y1;
		y_all_steps[step * 2 * gridDim.x + 1] = y2;
	}
}

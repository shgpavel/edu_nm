#include <math.h>

__global__ void rk4_cuda(float h, float A, float B,
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
		
		float k2_y1 = A * (y2 + 0.5 * h * k1_y1); 
		float k2_y2 = -B * (y1 + 0.5 * h * k1_y2);

		float k3_y1 = A * (y2 + 0.5 * h * k2_y1); 
		float k3_y2 = -B * (y1 + 0.5 * h * k2_y2);

		float k4_y1 = A * (y2 + h * k3_y1);
		float k4_y2 = -B * (y1 + h * k3_y2);

		/* float k_sum_1 = h * (A * y2
											 + 2 * A * (y2 + 0.5 * h * A * y2)
											 + 2 * A * (y2 + 0.5 * h * A * (y2 + 0.5 * h * A * y2))
												 + A * (y2 + h * A * (y2 + 0.5 * h * A * (y2 + 0.5 * h * A * y2))));
		*/
	
	  float k_sum_1 = h * (k1_y1 + 2 * k2_y1 + 2 * k3_y1 + k4_y1);
	  float k_sum_2 = h * (k1_y2 + 2 * k2_y2 + 2 * k3_y2 + k4_y2);
		
		y1 = y1 + k_sum_1 / 6;
		y2 = y2 + k_sum_2 / 6;
    
		y_all_steps[step * 2 * gridDim.x] = y1;
		y_all_steps[step * 2 * gridDim.x + 1] = y2;
	}
}

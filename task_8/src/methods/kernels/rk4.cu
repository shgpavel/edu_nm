#include <math.h>

__global__ void rk4_cuda(double h, double A, double B,
												 unsigned num_steps, double *y_initial,
												 double *y_values) {
    
	int system_idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	double *y_current = y_initial + 2 * system_idx;
	double *y_all_steps = y_values + system_idx;
	
	double y1 = y_current[0];
	double y2 = y_current[1];

	y_all_steps[0] = y1;
	y_all_steps[num_steps * 2 * gridDim.x] = y2;
	
	for (unsigned step = 1; step < num_steps; ++step) {
		double k1_y1 = A * y2;
		double k1_y2 = -B * y1;
		
		double k2_y1 = A * (y2 + 0.5 * h * k1_y1); 
		double k2_y2 = -B * (y1 + 0.5 * h * k1_y2);

		double k3_y1 = A * (y2 + 0.5 * h * k2_y1); 
		double k3_y2 = -B * (y1 + 0.5 * h * k2_y2);

		double k4_y1 = A * (y2 + h * k3_y1);
		double k4_y2 = -B * (y1 + h * k3_y2);

	  double k_sum_1 = h * (k1_y1 + 2 * k2_y1 + 2 * k3_y1 + k4_y1);
	  double k_sum_2 = h * (k1_y2 + 2 * k2_y2 + 2 * k3_y2 + k4_y2);
		
		y1 = y1 + k_sum_1 / 6;
		y2 = y2 + k_sum_2 / 6;
    
		y_all_steps[step * 2 * gridDim.x] = y1;
		y_all_steps[step * 2 * gridDim.x + 1] = y2;
	}
}

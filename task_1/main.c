/*
Copyright 2023 Pavel Shago

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include <stdio.h>
#include <math.h>

#include "funcs.h"

#define epsilon_sigma 1e-14

void z(double x, double epsilon_u,
       double epsilon_phi, double epsilon_psi, 
                FILE *csvout) {
    
    double u = sinus(power(x, 2) + 0.4, epsilon_u);
    double phi = square_rootPD(1 + power(x, 2), epsilon_phi) / (1 - x);
    double psi = sinhus(phi, epsilon_psi);
    
    double u_ = sin(pow(x, 2) + 0.4);
    double phi_ = sqrt(pow(x, 2) + 1) / (1 - x);
    double psi_ = sinh(phi_);

    fprintf(csvout, "%lf, %lf, %e, %lf, %e, %lf, %e, %lf, %e, %lf, "
            "%e, %lf, %e, %lf, %e, %lf, %lf\n",
            x, u, epsilon_u, u_, fabs(u_ - u), phi, epsilon_phi, phi_,
            fabs(phi_ - phi), psi, epsilon_psi, psi_, fabs(psi_ - psi),
            psi/u, 1e-6, psi_/u_, fabs(psi_/u_ - psi/u));
}

int main(void) {
    double lower, upper, step;
    double epsilon_u, epsilon_phi, epsilon_psi;
    FILE *csvout = fopen("output.csv", "w");

    printf( "Provide the input range as two floats:\n" );
    scanf( "%lf%lf", &lower, &upper );

    printf( "Provide step info:\n" );
    scanf( "%lf", &step );

    printf( "Provide epsilons for u func, phi, psi:\n" );

    scanf( "%lf%lf%lf", &epsilon_u, &epsilon_phi, &epsilon_psi );
    
    fprintf(csvout, "x, u(x), delta_u, u_(x), delta_u_, phi(x), delta_phi, phi_(x), "
                    "delta_phi_, psi(x), delta_psi, psi_(x), delta_psi_, z(x), "
                    "delta_z, z_(x), delta_z_\n");

    for (double i = lower; i < upper + epsilon_sigma; i += step) {
        z(i, epsilon_u, epsilon_phi, epsilon_psi, csvout);
    }

    printf( "Data saved in output.csv\n" );
    
    fclose(csvout);
    
    return 0;
}
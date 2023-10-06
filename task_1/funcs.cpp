#include <iostream>

double factorial(int n) {
    if (n <= 1) return 1.0;
    return n * factorial(n - 1);
}

double pow(double base, int exponent) {
    double result = 1.0;
    if (exponent > 0) {
        for (int i = 0; i < exponent; i++) {
            result *= base;
        }
    } else if (exponent < 0) {
        for (int i = 0; i > exponent; i--) {
            result /= base;
        }
    }
    return result;
}

double sin(double x, int n) {
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        result += (pow(-1, i) / factorial(2 * i + 1)) * pow(x, 2 * i + 1);
    }
    return result;
}

double sinh(double x, int n) {
    double result = 0.0;
    for (int i = 0; i < n; i++) {
      result += pow(x, 2 * i + 1) / factorial(2 * i + 1);
    }
    return result;
}

double square_rootPD(double num, double epsilon) {

    double guess = num;
    double prevguess = 0.0;

    while (std::abs(guess - prevguess) > epsilon) {
        prevguess = guess;
        guess = 0.5 * (guess + num / guess);
    }

    return guess;
}
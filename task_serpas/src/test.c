#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>
#include <math.h>

int
main()
{
	void    *handle = dlopen("./sca25lib.so", RTLD_LAZY);
	double  *(*funcTask1)(double **a, double *y, int N);
	long int (*funcTask2)(int **M, int N);
	double **(*funcTask3)(float **A, float **B, int N);

	funcTask1 = (double *(*)(double **, double *, int))dlsym(handle,
	                                                         "funcTask1");

	funcTask2 = (long int (*)(int **, int))dlsym(handle, "funcTask2");

	funcTask3 = (double **(*)(float **, float **, int))dlsym(handle,
	                                                         "funcTask3");

	/*
	        // 1
	        int n;
	        scanf("%d", &n);

	        double **a = (double **)malloc(sizeof(double *) * n);
	        for (int i = 0; i < n; ++i) {
	                a[i] = (double *)malloc(sizeof(double) * n);
	        }

	        for (int i = 0; i < n; ++i) {
	                for (int j = 0; j < n; ++j) {
	                        scanf("%lf", &a[i][j]);
	                }
	        }

	        double *y = (double *)malloc(sizeof(double) * n);

	        for (int i = 0; i < n; ++i) {
	                scanf("%lf", &y[i]);
	        }

	        double *res = funcTask1(a, y, n);
	        for (int i = 0; i < n; ++i) {
	                printf("%lg ", res[i]);
	        }
	        printf("\n");

	        double test_1[] = { 1, 2, 3 };
	        size_t flag     = 1;
	        for (int i = 0; i < n; ++i) {
	                if (fabs(test_1[i] - res[i]) > 1e-5)
	                        flag = 0;
	        }

	        flag == 1 ? printf("PASSED\n") : printf("FAILED\n");

	        free(res);

	        // 2-5
	        for (size_t k = 2; k < 6; ++k) {
	                scanf("%d", &n);

	                for (int i = 0; i < n; ++i) {
	                        for (int j = 0; j < n; ++j) {
	                                scanf("%lf", &a[i][j]);
	                        }
	                }

	                for (int i = 0; i < n; ++i) {
	                        scanf("%lf", &y[i]);
	                }

	                res = funcTask1(a, y, n);
	                for (int i = 0; i < n; ++i) {
	                        printf("%lg ", res[i]);
	                }
	                printf("\n");

	                double  test_2[] = { 1, 1, 1 };
	                double  test_3[] = { 1.26175, 1.23557, 1.21416 };
	                double  test_4[] = { 1.23704, 1.09637, 1.06815 };
	                double  test_5[] = { -0.174935, 0.608355, 0.765013 };
	                double *tests[]  = { test_2, test_3, test_4, test_5 };

	                flag             = 1;
	                for (int i = 0; i < n; ++i) {
	                        if (fabs(tests[k - 2][i] - res[i]) > 1e-5)
	                                flag = 0;
	                }

	                flag == 1 ? printf("PASSED\n") : printf("FAILED\n");

	                free(res);
	        }


	                for (int i = 0; i < n; ++i) {
	                  for (int j = 0; j < n; ++j) {
	                          printf("%lg ", a[i][j]);
	                  }
	                  printf("\n");
	                }

	                for (int i = 0; i < n; ++i) {
	                        printf("%lg ", y[i]);
	                }
	                printf("\n");


	        for (int i = 0; i < n; ++i) {
	                free(a[i]);
	        }
	        free(a);

	        free(y);
	*/

	// 2

	int n;
	scanf("%d", &n);

	int **a = (int **)malloc(sizeof(int *) * n);
	for (size_t i = 0; i < n; ++i) {
		a[i] = (int *)malloc(sizeof(int) * n);
	}

	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			scanf("%d", &a[i][j]);
		}
	}

	printf("%ld\n", funcTask2(a, n));

	for (size_t i = 0; i < n; ++i) {
		free(a[i]);
	}
	free(a);

	/*
	        // 3

	        int n;
	        scanf("%d", &n);

	        float **a = (float **)malloc(sizeof(float *) * n);
	        float **b = (float **)malloc(sizeof(float *) * n);
	        for (int i = 0; i < n; ++i) {
	                a[i] = (float *)malloc(sizeof(float) * n);
	                b[i] = (float *)malloc(sizeof(float) * n);
	        }

	        for (int i = 0; i < n; ++i) {
	                for (int j = 0; j < n; ++j) {
	                        scanf("%f", &a[i][j]);
	                }
	        }

	        for (int i = 0; i < n; ++i) {
	                for (int j = 0; j < n; ++j) {
	                        scanf("%f", &b[i][j]);
	                }
	        }


	        double **c = funcTask3(a, b, n);
	        for (size_t i = 0; i < n; ++i) {
	          for (size_t j = 0; j < n; ++j) {
	                  printf("%lg ", c[i][j]);
	          }
	          printf("\n");
	        }

	        for (size_t i = 0; i < n; ++i) {
	                free(c[i]);
	        }
	        free(c);
	*/

	dlclose(handle);
}

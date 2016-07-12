#include "matrix_util.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html */
void dgetrf_(int *M, int *N, double *A, int *LDA, int*IPIV, int *INFO);

/* http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f.html */
void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA, int *IPIV,
             double *B, int *LDB, int *INFO);

double f(double x, double y){
    return 0;
}

int main(int argc, char** argv) {
    double Re = 10.0; // レイノルズ数
    double U0 = 1.0;  // 初期速度
    double delta_t = 0.01;
    double x_width = 0.2;
    double y_width = 0.2;
    int max_count = 3;
    double **u;
    double **v;
    double **u_p;
    double **v_p;
    double **u_next;
    double **v_next;
    double **a;
    double **phi;
    double *b;
    int *ipiv;
    double **p;
    double h;

    int info;
    char trans = 'T';
    int nrhs = 1;

    int n, count;
    int i, j, k;

    if (argc < 5 + 1) {
        fprintf(stderr, "Usage : %s n Re U0 delta_t max_count\n", argv[0]);
    }
    n = atoi(argv[1]);
    Re = atof(argv[2]);
    U0 = atof(argv[3]);
    delta_t = atof(argv[4]);
    max_count = atoi(argv[5]);

    u = alloc_dmatrix(n+1, n+1);
    v = alloc_dmatrix(n+1, n+1);
    p = alloc_dmatrix(n+1, n+1);
    u_p = alloc_dmatrix(n+1, n+1);
    v_p = alloc_dmatrix(n+1, n+1);
    u_next = alloc_dmatrix(n+1, n+1);
    v_next = alloc_dmatrix(n+1, n+1);
    h = 1.0 / ((double)n);
    int m = (n+1) * (n+1);
    int x, y;
    a = alloc_dmatrix(m, m);
    phi = alloc_dmatrix(n+1, n+1);
    b = alloc_dvector(m);
    ipiv = alloc_ivector(m);

    // init
    for (i = 0; i <= n; ++i) {
        for (j = 0; j <= n; ++j) {
            v[i][j] = 0;
			v_p[i][j] = 0;

            if (i == 0 || i == n || j == 0 || j == n || j * h <= 0.25 || j * h >= 0.75 || i * h < 0.2) {
				u[i][j] = U0;
				u_p[i][j] = U0;
            } else if (0.5 - x_width / 2.0 <= i * h && i * h <= 0.5 + x_width / 2.0 && 
                       0.5 - y_width / 2.0 <= j * h && j * h <= 0.5 + y_width / 2.0 ) {
                    u_p[i][j] = 0.0;
                    v_p[i][j] = 0.0;
                    continue;
            } else {
				u[i][j] = 0;
				u_p[i][j] = 0;
			}

            phi[i][j] = 0;
			p[i][j] = 0;
        }
    }

    for (k = 0; k < m; ++k) {
        // a : 係数行列
        x = (k % (n+1) );
        y = (k / (n+1) );

		if (x == 0 || x == n || y == 0 || y == n) {
			a[k][k] = 1.0;
		} else {
        	a[k][k] = -4.0;
        	a[k][k - 1] = 1.0;
			a[k][k + 1] = 1.0;
			a[k][k - (n+1)] = 1.0;
			a[k][k + (n+1)] = 1.0;
			/*
        	if (x != 0) a[k][k - 1] = 1.0;
        	if (x != n) a[k][k + 1] = 1.0;
        	if (y != 0) a[k][k - (n+1)] = 1.0;
        	if (y != n) a[k][k + (n+1)] = 1.0;
			*/
		}
    }

    /* perform LU decomposition */
    ipiv = alloc_ivector(m);
    dgetrf_(&m, &m, &a[0][0], &m, &ipiv[0], &info);
    if (info != 0) {
        fprintf(stderr, "Error: LAPACK::dgetrf failed\n");
        exit(1);
    }

    // calculation
    count = 0;
    while (count < max_count) {
        // 中間流速u_p,v_p
        for (i = 0; i <=n; ++i) {
            for (j = 0; j <= n; ++j) {
                if (i == 0 || i == n || j == 0 || j == n ) {
                    u_p[i][j] = U0;
					v_p[i][j] = 0.0;
                    continue;
                } else if (0.5 - x_width / 2.0 <= i * h && i * h <= 0.5 + x_width / 2.0 && 
                           0.5 - y_width / 2.0 <= j * h && j * h <= 0.5 + y_width / 2.0 ) {
                    u_p[i][j] = 0.0;
                    v_p[i][j] = 0.0;
                    continue;
                }

                u_p[i][j] = -(p[i+1][j] - p[i-1][j]) / (2*h)
                            + u[i][j] * (- (u[i+1][j] - u[i-1][j]) / (2*h) - 4 / (Re * h * h))
                            + v[i][j] * (- (u[i][j+1] - u[i][j-1]) / (2*h))
                            + (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1]) / (Re * h * h);
                u_p[i][j] = delta_t * u_p[i][j] + u[i][j];
                
                v_p[i][j] = -(p[i][j+1] - p[i][j-1]) / (2*h)  
                            + v[i][j] * (- (v[i][j+1] - v[i][j-1]) / (2*h) - 4 / (Re * h * h))
                            + u[i][j] * (- (v[i+1][j] - v[i-1][j]) / (2*h))
                            + (v[i+1][j] + v[i-1][j] + v[i][j+1] + v[i][j-1]) / (Re * h * h);
                v_p[i][j] = delta_t * v_p[i][j] + v[i][j];
            }
        }

        // phi
        for (k = 0; k < m; ++k) {
            x = (k % (n+1) );
            y = (k / (n+1) );
			if (x == 0 || x == n || y == 0 || y == n){
				b[k] = 0.0;
			} else {
            	b[k] = u_p[x+1][y] - u_p[x-1][y] + v_p[x][y+1] - v_p[x][y-1];
            	b[k] = (1.0) * b[k] * h / delta_t / 2.0;
			}
        }
        /* solve equations */
        dgetrs_(&trans, &m, &nrhs, &a[0][0], &m, &ipiv[0], &b[0], &m, &info);
        if (info != 0) {
            fprintf(stderr, "Error: LAPACK::dgetrs failed\n");
            exit(1);
        }
        for (k = 0; k < m; ++k) {
            x = (k % (n + 1));
            y = (k / (n + 1));
            phi[x][y] = b[k];
        }

        for (i = 1; i <n; ++i) {
            for (j = 1; j < n; ++j) {
                if (i == 0 || i == n || j == 0 || j == n ) {
                    u[i][j] = U0;
                    continue;
                } else if (0.5 - x_width / 2.0 <= i * h && i * h <= 0.5 + x_width / 2.0 && 
                           0.5 - y_width / 2.0 <= j * h && j * h <= 0.5 + y_width / 2.0 ) {
                    u[i][j] = 0.0;
                    v[i][j] = 0.0;
                    continue;
                }


                u[i][j] = u_p[i][j] - delta_t * (phi[i+1][j] - phi[i-1][j]) / (2*h);
                v[i][j] = v_p[i][j] - delta_t * (phi[i][j+1] - phi[i][j-1]) / (2*h);

                p[i][j] = p[i][j] + phi[i][j];
            }
        }

        ++count;
    }

    printf("# x y u v \n");
    for (i = 0; i <= n; ++i) {
        for (j = 0; j <= n; ++j) {
            printf("%lf %lf %lf %lf\n", i * h, j * h, u[i][j], v[i][j]);
        }
        printf("\n");
    }


    free_dmatrix(u);
    free_dmatrix(v);
    free_dmatrix(p);
    free_dmatrix(u_p);
    free_dmatrix(v_p);
    free_dmatrix(u_next);
    free_dmatrix(v_next);
    free_dmatrix(a);
    free_dmatrix(phi);
    free_ivector(ipiv);
    free_dvector(b);    


    return 0;
}

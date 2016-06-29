#include "matrix_util.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char** argv) {
    double Re = 10.0; // レイノルズ数
    double U0 = 1.0;  // 初期速度
    double delta_t = 0.01;
    double x_width = 0.2;
    double y_width = 0.2;
    int max_count = 3;
    double **u;
    double **v;
    double **u_next;
    double **v_next;
    double h;
    
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


    u = alloc_dmatrix(n+2, n+2);
    v = alloc_dmatrix(n+2, n+2);
    u_next = alloc_dmatrix(n+2, n+2);
    v_next = alloc_dmatrix(n+2, n+2);
    h = 1.0 / ((double)n + 1.0);

    // init
    for (i = 0; i < n+2; ++i) {
        for (j = 0; j < n+2; ++j) {
            v[i][j] = 0;

            if (i == 1 || i == n || j == 1 || j == n ) u[i][j] = U0;
            else u[i][j] = 0;
        }
    }

    count = 0;
    while (count < max_count) {
        for (i = 1; i <=n; ++i) {
            for (j = 1; j <= n; ++j) {
                if (i == 1 || i == n || j == 1 || j == n ) {
                    u_next[i][j] = U0;
                    continue;
                } else if (0.5 - x_width / 2.0 <= i * h && i * h <= 0.5 + x_width / 2.0 && 
                           0.5 - y_width / 2.0 <= j * h && j * h <= 0.5 + y_width / 2.0 ) {
                    u_next[i][j] = 0.0;
                    v_next[i][j] = 0.0;
                    continue;
                }
                u_next[i][j] =   u[i][j] * (- (u[i+1][j] - u[i-1][j]) / (2*h) - 4 / (Re * h * h))
                               + v[i][j] * (- (u[i][j+1] - u[i][j-1]) / (2*h))
                               + (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1]) / (Re * h * h);
                u_next[i][j] = delta_t * u_next[i][j] + u[i][j];
                
                v_next[i][j] =   v[i][j] * (- (v[i][j+1] - v[i][j+1]) / (2*h) - 4 / (Re * h * h))
                               + u[i][j] * (- (u[i+1][j] - u[i-1][j]) / (2*h))
                               + (v[i+1][j] + v[i-1][j] + v[i][j+1] + v[i][j-1]) / (Re * h * h);
                v_next[i][j] = delta_t * v_next[i][j] + v[i][j];
            
            }
        }

        for (i = 1; i <=n; ++i) {
            for (j = 1; j <= n; ++j) {
                u[i][j] = u_next[i][j];
                v[i][j] = v_next[i][j];
            }
        }

        ++count;
    }

    printf("# x y u v \n");
    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j) {
            printf("%lf %lf %lf %lf\n", i * h, j * h, u[i][j], v[i][j]);
        }
        printf("\n");
    }


    free_dmatrix(u);
    free_dmatrix(v);
    free_dmatrix(u_next);
    free_dmatrix(v_next);




    return 0;
}

//Header file for OMP
#ifndef _OMP_H_
#define _OMP_H_

//defining the dimensions of operating matrices
#define P 8
#define L 299
#define M 8
#define N 14
#define K 8

void INIT_A1_hw(float A[P*N], float A_i[P*N],float A1[P*N], float A1_i[P*N]);
float norm_hw(float matrix[P*N], int j);
void update_Z_hw(float A1[P*N], float A1_i[P*N],float R[P*L], float R_i[P*L],float Z[N]);
void complexMatrixMult_hw(float matrix_1[K*P], float matrix_1_i[K*P], float matrix_2[P*L], float matrix_2_i[P*L], int p, int l, int m, int n, float matrix_product[K*L], float matrix_product_i[K*L]);
#endif

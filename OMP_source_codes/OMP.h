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
void update_As_hw(float A1[P*N],float As[P*K],int k, int maxind);
void update_Z_hw(float A1[P*N], float A1_i[P*N],float R[P*L], float R_i[P*L],float Z[N]);
void complexMatrixMult_hw(float matrix_1[K*P], float matrix_1_i[K*P], float matrix_2[P*L], float matrix_2_i[P*L], int p, int l, int m, int n, float matrix_product[K*L], float matrix_product_i[K*L]);
void get_Ax_hw(float As[P*N], float As_i[P*N],float X[K*L], float X_i[K*L],int k, float Ax[P*L], float Ax_i[P*L]);
void update_R_hw(float R[P*L], float R_i[P*L],float Ax[P*L],float Ax_i[P*L]);

////////////////////////////////////////////////////////////////////////////////////////////

void MatrixMult_(float matrix_1[K*K], float matrix_2[K*K], int p, int l, int m, int n, float matrix_product[K*K]);
bool lup_(float A[K*K], float Lo[K*K], float U[K*K], float Pe[K], int DIM);
void Lower_inv_(float Lo[K*K], float L_inv[K*K], int DIM);
void Upper_inv_(float U[K*K], float U_inv[K*K], int DIM);
void matrix_mult_(float U_inv[K*K], float L_inv[K*K], float A_inv[K*K], int DIM);
void final_perm_(float UL_inv[K*K], float Pe[K*K], float A_inv[K*K], int DIM);
bool inverse_(float A[K*K], float A_inv[K*K], int DIM);
void Add_(float mat1[K*K], float mat2[K*K], int p, int l, float mat12[K*K]);
void neg_(float neg_inv[K*K], int m,int n, float inv[K*K]);
//finds the complex inverse of matrix
void complex_inverse_(float A[K*K], float Ai[K*K], int n, float inv[K*K], float inv_i[K*K]);
//finds the pseudo inverse of complex matrix
void pseudo_(float matrix[K*K],float matrix_i[K*K], int p, int k, float pinv[K*K], float pinv_i[K*K]);


#endif

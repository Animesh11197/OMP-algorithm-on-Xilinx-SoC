// contains all the potential HW Functions:
#include "OMP.h"
#include <iostream>
#include <cstdio>
#include<stdio.h>
#include<stdlib.h>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <inttypes.h>
#include <sds_lib.h>
#include "sds_lib.h"
#include "ap_int.h"
#include <ap_fixed.h>
using namespace std;
#include <math.h>


float norm_hw(float matrix[P*N], int j)
{
	float sum = 0;
	for (int i = 0; i < P; i++){
#pragma HLS PIPELINE

		sum = sum + (matrix[i*N + j] * matrix[i*N + j]);}
	return sqrt(sum);
}

void INIT_A1_hw(float A[P*N], float A_i[P*N],float A1[P*N], float A1_i[P*N])
{
	//absolute matrix of A i.e. A+Ai
	float abs_A[P*N];
	for (int i = 0; i < P*N; i++)
	{
		#pragma HLS PIPELINE

		abs_A[i] = sqrt(A[i] * A[i] + A_i[i] * A_i[i]);
	}
	//display(abs_A,P,N);
	for (int j = 0; j < N; j++)
	{
		#pragma HLS PIPELINE

		float b = norm_hw(abs_A, j);
		//printf("\nNorm(A,j)= %f\t %d\n", b, j);
		for (int i = 0; i < P; i++)
		{
			A1[i*N + j] = A[i*N + j] / b;
			A1_i[i*N + j] = A_i[i*N + j] / b;
		}
	}
}

// initializing new rows in As at every iteration of k
void update_As_hw(float A1[P*N],float As[P*K],int k, int maxind)
{

	for (int i = 0; i < P; i++)
	{
#pragma HLS UNROLL

		As[i*K + k] = A1[i*N + maxind];

	}
}
void update_Z_hw(float A1[P*N], float A1_i[P*N],float R[P*L], float R_i[P*L],float Z[N])
{
	float term1,term2,term3,term4;
	float product,product_i, sum,sum_i;
	int index1,index2;
	// Loop for updating Z
	for (int j = 0; j < N; j++)
	{
//#pragma HLS PIPELINE
		sum = 0;
		for (int i = 0; i < L; i++)
		{
			product = 0; product_i = 0;
			for (int z = 0; z < P; z++)
			{
#pragma HLS PIPELINE
				index1=z*N + j;
				index2=z*L + i;
				term1=(A1[index1] * R[index2]);
				term2=(-A1_i[index1] * R_i[index2]);
				term3=(A1[index1] * R_i[index2]);
				term4=(-A1_i[index1] * R[index2]);
				product += (term1-term2);
				product_i += (term3 + term4);
			}
			sum += (product*product)+ (product_i*product_i);
		}
		//cout << endl;
		Z[j] = sqrt(sum);
	}
}

//complex matrix multiplication function
void complexMatrixMult_hw(float matrix_1[K*P], float matrix_1_i[K*P], float matrix_2[P*L], float matrix_2_i[P*L], int p, int l, int m, int n, float matrix_product[K*L], float matrix_product_i[K*L])
{
	for (int i = 0; i < p; i++) {
#pragma HLS loop_tripcount max=8
//#pragma HLS PIPELINE
		for (int j = 0; j < n; j++) {             // not j<M
#pragma HLS loop_tripcount max=299
			matrix_product[i * n + j] = 0;
			matrix_product_i[i * n + j] = 0;
			for (int k = 0; k < l; k++) {
#pragma HLS loop_tripcount max=8
#pragma HLS PIPELINE
				matrix_product[i * n + j] += (matrix_1[i * l + k] * matrix_2[k * n + j]) - (matrix_1_i[i * l + k] * matrix_2_i[k * n + j]);
				matrix_product_i[i * n + j] += (matrix_1[i * l + k] * matrix_2_i[k * n + j]) + (matrix_1_i[i * l + k] * matrix_2[k * n + j]);;
			}
		}
	}
	return;

}

/*void complexMatrixMult_hw(float matrix_1[K*P], float matrix_1_i[K*P], float matrix_2[P*L], float matrix_2_i[P*L], int p, int l, int m, int n, float matrix_product[K*L], float matrix_product_i[K*L])
{
	float local_matrix_1[K][P],local_matrix_2[P][L],local_matrix_1_i[K][P],local_matrix_2_i[P][L];

#pragma HLS ARRAY_PARTITION variable= local_matrix_1 complete dim=2
#pragma HLS ARRAY_PARTITION variable= local_matrix_1_i complete dim=2
#pragma HLS ARRAY_PARTITION variable= local_matrix_2 complete dim=1
#pragma HLS ARRAY_PARTITION variable= local_matrix_2_i complete dim=1

			for(int i=0;i<K;i++)
			{
#pragma HLS PIPELINE
				for(int j=0;j<P;j++){

					local_matrix_1[i][j]=matrix_1[i*P+j];
					local_matrix_1_i[i][j]=matrix_1_i[i*P+j];}
			}
			for(int i=0;i<P;i++)
				{
#pragma HLS PIPELINE
				for(int j=0;j<L;j++)
					{

					local_matrix_2[i][j]=matrix_2[i*L+j];
					local_matrix_2_i[i][j]=matrix_2_i[i*L+j];			}}



	for (int i = 0; i < K; i++) {
#pragma HLS PIPELINE
		for (int j = 0; j < L; j++) {
//#pragma HLS PIPELINE
			float result=0;
			float result2=0;

			for (int k = 0; k < P; k++) {
//#pragma HLS PIPELINE
				float term1= (local_matrix_1[i][k]);
				float term2= (local_matrix_2[k][j]);
				float term3= (local_matrix_1_i[i][k]);
				float term4= (local_matrix_2_i[k][j]);
				result +=  (term1*term2)-(term3*term4);
				result2 +=  (term1*term4)+(term3*term2);

				//matrix_product_i[i * L + j] += (local_matrix_1[i][k] * local_matrix_2_i[k][j]) + (local_matrix_1_i[i][k] * local_matrix_2[k][j]);
			}

			matrix_product[i * L + j]=result;
			matrix_product_i[i * L + j]=result2;


		}
	}


	return;

}*/

//get Ax=As*X
void get_Ax_hw(float As[P*N], float As_i[P*N],float X[K*L], float X_i[K*L],int k, float Ax[P*L], float Ax_i[P*L])
{
	//MatrixMult(As, X, P, k + 1, k + 1, L, Ax);
	for (int i = 0; i < P; i++)
	{
		for (int j = 0; j < L; j++)
		{
			Ax[i * L + j] = 0; Ax_i[i * L + j] = 0;
			for (int z = 0; z < (k + 1); z++)
			{
#pragma HLS loop_tripcount max=8
#pragma HLS PIPELINE
				Ax[i * L + j] += (As[i * K + z] * X[z * L + j]) - (As_i[i * K + z] * X_i[z * L + j]);
				Ax_i[i * L + j] += (As_i[i * K + z] * X[z * L + j]) + (As[i * K + z] * X_i[z * L + j]);
			}
		}
	}

}

//UPDATING RESDUE
void update_R_hw(float R[P*L], float R_i[P*L],float Ax[P*L],float Ax_i[P*L])
{
	int index;
	//UPDATING RESIDUE MATRIX
	for (int i = 0; i < P; i++)
	{
		for (int j = 0; j < L; j++)
		{
#pragma HLS PIPELINE

			index=i*L + j;
			R[index] = R[index] - Ax[index];
			R_i[index] = R_i[index] - Ax_i[index];
		}
	}
}
//Add complex pseudo inverse hw function


//function to multiply Matrices
void MatrixMult_(float matrix_1[K*K], float matrix_2[K*K], int p, int l, int m, int n, float matrix_product[K*K]) {
	for (int i = 0; i < p; i++) {
#pragma HLS loop_tripcount max=8
		for (int j = 0; j < n; j++) {             // not j<M
#pragma HLS loop_tripcount max=8
			matrix_product[i * n + j] = 0;
			for (int k = 0; k < l; k++) {
#pragma HLS loop_tripcount max=8
#pragma HLS PIPELINE
				matrix_product[i * n + j] += matrix_1[i * l + k] * matrix_2[k * n + j];
			}
		}
	}
	return;

}


bool lup_(float A[K*K], float Lo[K*K], float U[K*K], float Pe[K], int DIM)
{
PermutMat_Initialize: for (int i = 0; i < DIM; i++)
{
#pragma HLS loop_tripcount max=8
	//#pragma HLS PIPELINE
	Pe[i] = i;
}

					  int i, j, k;
					  bool singular = 0;		/*-1 if matrix is singular --> inverse doesn't exist.*/

					  /*LUP Decomposition is finding L,U,P such that
					   * PA=LU, here
					   * P-Row Permutation Matrix
					   * A-Input Matrix
					   * L-Lower Triangular Matrix
					   * U-Upper Triangular Matrix
					   * Then the inverse for A is-
					   * inv(A)=inv(U)*inv(L)*P
					   * */
				  lup_label0: for (i = 0; i < DIM; i++)
				  {
#pragma HLS loop_tripcount max=8
					  //#pragma HLS PIPELINE
											/* pos-index for maximum value in a column.*/
					  int pos = i;

					  /*stores the maximum value, initially taken as follows*/
					  float max = A[i*DIM + i];

					  /*pivoting
					   * Find the maximum value in the column below the diagonal element (for the column)
					   * Swap the row with the one with max value.
					   * */
				  find_max: for (k = i + 1; k < DIM; k++)
				  {
#pragma HLS loop_tripcount max=8
					  //#pragma HLS PIPELINE
					  //#pragma HLS LOOP_TRIPCOUNT min=1 max=3
					  float tmp = A[k*DIM + i];
					  if (tmp > A[i*DIM + i] & tmp > max)
					  {
						  pos = k;
						  max = tmp;
					  }
				  }
							/*check if the matrix is singular
							 * */
							if (A[pos*DIM + i] == 0.0)
							{
								singular = 1; 		//matrix is singular
								return singular;
							}
							else
							{
								if (pos != i)
								{
								swap_row: for (k = 0; k < DIM; k++)
								{
									#pragma HLS loop_tripcount max=8

									float tmp = A[pos*DIM + k];
									A[pos*DIM + k] = A[i*DIM + k];
									A[i*DIM + k] = tmp;
								}

										  /*update the permutation matrix for the swap*/
										  int ind1, ind2;
										  for (int i1 = 0; i1 < DIM; i1++)
										  {
											  #pragma HLS loop_tripcount max=8
											  if (Pe[i1] == pos)
												  ind1 = i1;
											  if (Pe[i1] == i)
												  ind2 = i1;
										  }
										  float temp = Pe[ind2];
										  Pe[ind2] = Pe[ind1];
										  Pe[ind1] = temp;
								}
							}

							/*extract the L and U Matrices
							 * */
						lup_label1: for (k = i + 1; k < DIM; k++)
						{
#pragma HLS loop_tripcount max=8
							A[k*DIM + i] = A[k*DIM + i] / A[i*DIM + i];
						lup_label2: for (j = i + 1; j < DIM; j++)
						{
#pragma HLS loop_tripcount max=8
							A[k*DIM + j] = A[k*DIM + j] - A[i*DIM + j] * A[k*DIM + i];
						}
						}
				  }
							  //L matrix: Lower half of A matrix and diagonal elements is 1.
						  Assign_L0: for (int i = 0; i < DIM; i++)
						  {
#pragma HLS loop_tripcount max=8
						  Assign_L1: for (int j = i; j < DIM; j++)
						  {
#pragma HLS loop_tripcount max=8
							  if (i == j)
								  Lo[j*DIM + i] = 1;
							  else
								  Lo[j*DIM + i] = A[j*DIM + i];
						  }
						  }

									 //U matrix: Upper half of A matrix.
								 Assign_U0: for (int i = 0; i < DIM; i++)
								 {
#pragma HLS loop_tripcount max=8
									Assign_U1: for (int j = i; j < DIM; j++)
									{
#pragma HLS loop_tripcount max=8
										U[i*DIM + j] = A[i*DIM + j];
									}
								 }

											Pe[DIM] = (singular == 1) ? -1 : 0;
											return 0;
}

void Lower_inv_(float Lo[K*K], float L_inv[K*K], int DIM)
{
	/* To calculate the inverse of L matrix:
	 * We have,
	 * for i==j, L_inv(i,j)=1/L(i,j)
	 * 						  k=i
	 * for i>j,  L_inv(i,j)=-Summation{L(i,k)*L_inv(k,j)}
	 * 						  k=j
	 * */

linv_label0: for (int i = 0; i < DIM; i++)
{
#pragma HLS loop_tripcount max=8
linv_label1: for (int j = 0; j < DIM; j++)
{
#pragma HLS loop_tripcount max=8
	if (i < j)
		L_inv[i*DIM + j] = 0;
	else if (i == j)
		L_inv[i*DIM + j] = 1 / Lo[i*DIM + j];
	else
	{
		float sum = 0.0f;
	linv_label2: for (int k = j; k < i; k++)
					{
					#pragma HLS LOOP_TRIPCOUNT min=1 max=8

					//#pragma HLS PIPELINE
					sum = sum + Lo[i*DIM + k] * L_inv[k*DIM + j];
					}

				 L_inv[i*DIM + j] = -sum;
	}
}
}
}

void Upper_inv_(float U[K*K], float U_inv[K*K], int DIM)
{
	/* Inverse Of Upper Triangular Matrix
	 * We have,
	 * for i==j, U(i,j)=1/U(i,j)
	 * 						  k=i
	 * for i>j   U_inv(j,i)=[-Summation{L(i,k)*L_inv(k,j)} ]/U[i][i]
	 * 						  k=j
	 **/
uinv_label10: for (int i = 0; i < DIM; i++)
{
#pragma HLS loop_tripcount max=8
	//#pragma HLS PIPELINE
uinv_label11: for (int j = 0; j < DIM; j++)
{
#pragma HLS loop_tripcount max=8
	U_inv[i*DIM + j] = 0;
}
}

		  uinv_label20: for (int i = 0; i < DIM; i++)
		  {
#pragma HLS loop_tripcount max=8
			  //#pragma HLS PIPELINE
			  U_inv[i*DIM + i] = (float)1 / U[i*DIM + i];
		  }

					univ_label30: for (int i = 0; i < DIM; i++)
					{
#pragma HLS loop_tripcount max=8
					univ_label31: for (int j = 0; j < i; j++)
					{
#pragma HLS loop_tripcount max=8
						//#pragma HLS PIPELINE
						//#pragma HLS LOOP_TRIPCOUNT min=1 max=3
						float sum = 0.0f;
					univ_label32: for (int k = j; k < i; k++)
					{
#pragma HLS loop_tripcount max=8
						//#pragma HLS PIPELINE
						//#pragma HLS LOOP_TRIPCOUNT min=1 max=3
						sum += (U[k*DIM + i] * U_inv[j*DIM + k]);
					}
								  U_inv[j*DIM + i] = -sum / U[i*DIM + i];
					}
					}
}

void matrix_mult_(float U_inv[K*K], float L_inv[K*K], float A_inv[K*K], int DIM)
{
	float sumTemp = 0;
MM_L1: for (int i = 0; i < DIM; i++)
{
#pragma HLS loop_tripcount max=8
MM_L2: for (int j = 0; j < DIM; j++)
{
#pragma HLS loop_tripcount max=8
	//#pragma HLS PIPELINE
	float sumFinal = 0;
MM_L3: for (int k = 0; k < DIM; k++)
{
#pragma HLS loop_tripcount max=8
	sumTemp = U_inv[i*DIM + k] * L_inv[k*DIM + j];
	sumFinal += sumTemp;
}
	   A_inv[i*DIM + j] = sumFinal;
}
}
}

void final_perm_(float UL_inv[K*K], float Pe[K*K], float A_inv[K*K], int DIM)
{
		/*Multiplication with P is just the row permuation*/
L1: for (int i = 0; i < DIM; i++)
{
#pragma HLS loop_tripcount max=8
L2: for (int j = 0; j < DIM; j++)
{
#pragma HLS loop_tripcount max=8
	A_inv[j*DIM + i] = UL_inv[j*DIM + (int)Pe[i]];
}
}
}

bool inverse_(float A[K*K], float A_inv[K*K], int DIM)
{
	float Ac [K*K];
	for (int i = 0; i < DIM*DIM; i++)
	{
#pragma HLS loop_tripcount max=8*8
		Ac[i] = A[i];
	}
	float Lo [K*K];
	float U [K*K];
	float L_inv [K*K];
	float U_inv [K*K];
	float UL_inv [K*K];
	float Pe [K];

	// calculate the P,A,L,U such that PA=LU
	bool singular = lup_(Ac, Lo, U, Pe, DIM);

	if (singular)
	{
		cout << "SINGULAR\n";
		for (int i = 0; i < DIM; i++)
		{
#pragma HLS loop_tripcount max=8
			//#pragma HLS PIPELINE
			for (int j = 0; j < DIM; j++)
#pragma HLS loop_tripcount max=8
				A_inv[i*DIM + j] = 0;
		}
		return singular;
	}
	else
	{
		//cout << "///////////////////////////n";
		Lower_inv_(Lo, L_inv, DIM);
		Upper_inv_(U, U_inv, DIM);
		matrix_mult_(U_inv, L_inv, UL_inv, DIM);
		final_perm_(UL_inv, Pe, A_inv, DIM);


		/*printf("\n\nPrinting L_inverse Matrix\n");
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				printf("%f ", L_inv[i*DIM+j]);
			}
			printf("\n");
		}
		printf("\n\nPrinting U_inverse Matrix\n");
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				printf("%f ", U_inv[i*DIM+j]);
			}
			printf("\n");
		}
		printf("\n\nPrinting LU_inverse Matrix\n");
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				printf("%f ", UL_inv[i*DIM+j]);
			}
			printf("\n");
		}*/
		/*printf("\n\nPrinting A_inverse Matrix\n");
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				printf("%f ", A_inv[i*DIM + j]);
			}
			printf("\n");
		}*/

	}

	return 0;
}

//complex matrix multiplication function
/*void complexMatrixMult(float *matrix_1, float *matrix_1_i, float *matrix_2, float *matrix_2_i, int p, int l, int m, int n, float *matrix_product, float *matrix_product_i) {
	for (int i = 0; i < p; i++) {
		for (int j = 0; j < n; j++) {             // not j<M
			matrix_product[i * n + j] = 0;
			matrix_product_i[i * n + j] = 0;
			for (int k = 0; k < l; k++) {
				matrix_product[i * n + j] += (matrix_1[i * l + k] * matrix_2[k * n + j]) - (matrix_1_i[i * l + k] * matrix_2_i[k * n + j]);
				matrix_product_i[i * n + j] += (matrix_1[i * l + k] * matrix_2_i[k * n + j]) + (matrix_1_i[i * l + k] * matrix_2[k * n + j]);;
			}
		}
	}
	return;

}*/

// Adds matrix
void Add_(float mat1[K*K], float mat2[K*K], int p, int l, float mat12[K*K])
{
/*	for (int i = 0; i < p*l; i++)
	{
		mat12[i] = mat2[i] + mat1[i];
		//cout<<mat12[i]<<"  ";
	}*/

	for(int i=0;i<p;i++)
	{
#pragma HLS loop_tripcount max=8
		for(int j=0;j<l;j++)
		{
#pragma HLS loop_tripcount max=8
			mat12[i*l+j] = mat2[i*l+j] + mat1[i*l+j];
		}
	}
}

// Neg matrix function
void neg_(float neg_inv[K*K], int m,int n, float inv[K*K])
{
/*	for (int i = 0; i < m*n; i++)
	{
		inv[i] = -1 * neg_inv[i];
	}*/

	for(int i=0;i<m;i++)
	{
#pragma HLS loop_tripcount max=8
		for(int j=0;j<n;j++)
		{
#pragma HLS loop_tripcount max=8
			inv[i*n+j] = -1 * neg_inv[i*n+j];
		}
	}


}

//finds the complex inverse of matrix
void complex_inverse_(float A[K*K], float Ai[K*K], int n, float inv[K*K], float inv_i[K*K])
{
	float neg_inv_i [K*K];
	float A_i [K*K];
	float BA_i_B [K*K];
	float A_BA_i_B [K*K];
	float BA_i  [K*K];

	inverse_(A, A_i, n);
		//display(A_i, n, n);
	//display(A,n,n);
	MatrixMult_(A_i, Ai, n, n, n, n, BA_i);

	MatrixMult_(Ai, BA_i, n, n, n, n, BA_i_B);
	Add_(A, BA_i_B, n, n, A_BA_i_B);
	inverse_(A_BA_i_B, inv, n);
	//display(A_BA_i_B,n,n);
	//display(inv, n, n);
	MatrixMult_(BA_i, inv, n, n, n, n, neg_inv_i);
	neg_(neg_inv_i,n, n, inv_i);
	//display(inv_i, n, n);
}

//finds the pseudo inverse of complex matrix
void pseudo_(float matrix[8*8],float matrix_i[8*8], int p, int k, float pinv[8*8], float pinv_i[8*8])
{


	float matrix_mult[K*K]; 	float matrix_mult_i [K*K];

	for (int i = 0; i < k; i++)
	{
#pragma HLS loop_tripcount max=8
		for (int j = 0; j < k; j++)
		{
#pragma HLS loop_tripcount max=8
			matrix_mult[i * k + j] = 0;			matrix_mult_i[i * k + j] = 0;

			for (int z = 0; z < p; z++)
			{
#pragma HLS loop_tripcount max=8
#pragma HLS PIPELINE
				matrix_mult[i * k + j] += (matrix[z * K + i] * matrix[z * K + j])-(-matrix_i[z * K + i] * matrix_i[z * K + j]);
				matrix_mult_i[i * k + j] += (-matrix_i[z * K + i] * matrix[z * K + j]) + (matrix[z * K + i] * matrix_i[z * K + j]);

			}
			//cout << matrix_mult[i * K + j] << "\t" << matrix_mult_i[i * K + j]<<"---\n";
		}
	}
	float inv [K*K]; // To store inverse
	float inv_i [K*K]; // To store inverse
	//cout << "\n As'*As\n";
	//display(matrix_mult, k, k); display(matrix_mult_i, k, k);


	complex_inverse_(matrix_mult,matrix_mult_i,k,inv,inv_i);
	//cout << "\nINVERSE of square mat"<<k<<""<< ":\n";
	//display(inv, k, k); display(inv_i,k,k);


	//MatrixMult(inv, t_matrix, k, k, k, p, pinv);
	for (int i = 0; i < k; i++)
	{
#pragma HLS loop_tripcount max=8
		for (int j = 0; j < p; j++)
		{
#pragma HLS loop_tripcount max=8
			pinv[i * p + j] = 0;	pinv_i[i * p + j] = 0;
			for (int z = 0; z < k; z++)
			{
#pragma HLS loop_tripcount max=8
				pinv[i * p + j] += (inv[i * k + z] * matrix[j * p + z])+(inv_i[i * k + z] * matrix_i[j * p + z]);
				pinv_i[i * p + j] += ((inv_i[i * k + z] * matrix[j * p + z]) - (inv[i * k + z] * matrix_i[j * p + z]));
			}
		}
	}
	//cout << "\nThe Monroe-penrose inverse is :\n";
	//display(pinv, k, p);	display(pinv_i, k, p);



}

// contains all the potential HW Functions:
#include "OMP.h"
#include <math.h>
#include <ap_fixed.h>

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

void update_Z_hw(float A1[P*N], float A1_i[P*N],float R[P*L], float R_i[P*L],float Z[N])
{
	float product,product_i, sum,sum_i;

	// Loop for updating Z
	for (int j = 0; j < N; j++)
	{
#pragma HLS PIPELINE
		sum = 0;
		for (int i = 0; i < L; i++)
		{
			product = 0; product_i = 0;
			for (int z = 0; z < P; z++)
			{
				product += ((A1[z*N + j] * R[z*L + i])-(-A1_i[z*N + j] * R_i[z*L + i]));
				product_i += ((A1[z*N + j] * R_i[z*L + i]) + (-A1_i[z*N + j] * R[z*L + i]));
			}
			sum += (product*product)+ (product_i*product_i);
		}
		//cout << endl;
		Z[j] = sqrt(sum);
	}
}

void complexMatrixMult_hw(float matrix_1[K*P], float matrix_1_i[K*P], float matrix_2[P*L], float matrix_2_i[P*L], int p, int l, int m, int n, float matrix_product[K*L], float matrix_product_i[K*L])
{
	ap_fixed<47,40> local_matrix_1[K][P],local_matrix_2[P][L],local_matrix_1_i[K][P],local_matrix_2_i[P][L];

#pragma HLS ARRAY_PARTITION variable= local_matrix_1 complete dim=2
#pragma HLS ARRAY_PARTITION variable= local_matrix_1_i complete dim=2
#pragma HLS ARRAY_PARTITION variable= local_matrix_2 complete dim=1
#pragma HLS ARRAY_PARTITION variable= local_matrix_2_i complete dim=1

			for(int i=0;i<K;i++)
			{
				for(int j=0;j<P;j++){
#pragma HLS PIPELINE
					local_matrix_1[i][j]=matrix_1[i*P+j];
					local_matrix_1_i[i][j]=matrix_1_i[i*P+j];}
			}
			for(int i=0;i<P;i++)
				{
				for(int j=0;j<L;j++)
					{
#pragma HLS PIPELINE
					local_matrix_2[i][j]=matrix_2[i*L+j];
					local_matrix_2_i[i][j]=matrix_2_i[i*L+j];			}}



	for (int i = 0; i < K; i++) {
//#pragma HLS PIPELINE
		for (int j = 0; j < L; j++) {
#pragma HLS PIPELINE
			ap_fixed<47,40> result=0;
			ap_fixed<47,40> result2=0;

			for (int k = 0; k < P; k++) {
//#pragma HLS PIPELINE
				ap_fixed<47,40> term1= (local_matrix_1[i][k]);
				ap_fixed<47,40> term2= (local_matrix_2[k][j]);
				ap_fixed<47,40> term3= (local_matrix_1_i[i][k]);
				ap_fixed<47,40> term4= (local_matrix_2_i[k][j]);
				result +=  (term1*term2)-(term3*term4);
				result2 +=  (term1*term4)+(term3*term2);

				//matrix_product_i[i * L + j] += (local_matrix_1[i][k] * local_matrix_2_i[k][j]) + (local_matrix_1_i[i][k] * local_matrix_2[k][j]);
			}

			matrix_product[i * L + j]=result;
			matrix_product_i[i * L + j]=result2;


		}
	}


	return;

}

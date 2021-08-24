#include <iostream>
//#include <cstdlib>
#include <time.h>
using namespace std;
//#include "Header.h"

// C++ program to find Moore-Penrose inverse  matrix
#define P 8
#define L 299
#define M 8
#define N 14
#define K 8


// function to find row-wise normlisation of an array [MATLAB-->norm()]
float norm(float matrix[P*N], int j)
{
	float sum = 0;
	for (int i = 0; i < P; i++)
		sum = sum + (matrix[i*N + j] * matrix[i*N + j]);
	return sqrt(sum);
}

void Trans_2D_1D(float matrix_2D[P][N], float *matrix) {

	for (int i = 0; i < P; i++) {
		for (int j = 0; j < N; j++) {
			matrix[i * N + j] = matrix_2D[i][j];
		}
		//cout << endl;
	}

	return;

}




void Transpose(float *matrix, int p, int l, float *t_matrix) {
	for (int i = 0; i < p; i++) {
		for (int j = 0; j < l; j++)
		{
			t_matrix[j * p + i] = matrix[i * l + j];
		}
		//cout << endl;
	}
	return;
}


void MatrixMult(float *matrix_1, float *matrix_2, int p, int l, int m, int n, float *matrix_product) {
	for (int i = 0; i < p; i++) {
		for (int j = 0; j < n; j++) {             // not j<M
			matrix_product[i * n + j] = 0;
			for (int k = 0; k < l; k++) {
				matrix_product[i * n + j] += matrix_1[i * l + k] * matrix_2[k * n + j];
			}
		}
	}
	return;

}


bool lup(float *A, float *Lo, float *U, float *Pe, int DIM)
{

	//#pragma HLS INLINE
		/*Permutation matrix- Matrix obtained from the row permutation of the identity matrix.
		 * Initially, it will be taken as identity matrix. P initially assigned as [0 1 2]-implies the position of rows of identity matrix in P,
		 * as we swap rows in A matrix,we will update the P matrix simultaneously.\
		 * Suppose row0 swaped with row1 then P will be updated as [1 0 2] and the corresponding permutation matrix is 	0 1 0
		 * 																												1 0 0
		 * 																												0 0 1
		 * */
		 //Initialize Row Permutation Matrix
PermutMat_Initialize: for (int i = 0; i < DIM; i++)
{
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

									float tmp = A[pos*DIM + k];
									A[pos*DIM + k] = A[i*DIM + k];
									A[i*DIM + k] = tmp;
								}

										  /*update the permutation matrix for the swap*/
										  int ind1, ind2;
										  for (int i1 = 0; i1 < DIM; i1++)
										  {
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
							//#pragma HLS LOOP_TRIPCOUNT min=1 max=3
							A[k*DIM + i] = A[k*DIM + i] / A[i*DIM + i];
						lup_label2: for (j = i + 1; j < DIM; j++)
						{
							//#pragma HLS PIPELINE
							//#pragma HLS LOOP_TRIPCOUNT min=1 max=3
							A[k*DIM + j] = A[k*DIM + j] - A[i*DIM + j] * A[k*DIM + i];
						}
						}
				  }
							  //L matrix: Lower half of A matrix and diagonal elements is 1.
						  Assign_L0: for (int i = 0; i < DIM; i++)
						  {
						  Assign_L1: for (int j = i; j < DIM; j++)
						  {
							  //#pragma HLS PIPELINE
							  //#pragma HLS LOOP_TRIPCOUNT min=1 max=3
							  if (i == j)
								  Lo[j*DIM + i] = 1;
							  else
								  Lo[j*DIM + i] = A[j*DIM + i];
						  }
						  }

									 //U matrix: Upper half of A matrix.
								 Assign_U0: for (int i = 0; i < DIM; i++)
								 {
								 Assign_U1: for (int j = i; j < DIM; j++)
								 {
									 //#pragma HLS PIPELINE
									 //#pragma HLS LOOP_TRIPCOUNT min=1 max=3
									 U[i*DIM + j] = A[i*DIM + j];
								 }
								 }

											Pe[DIM] = (singular == 1) ? -1 : 0;
											return 0;
}

void Lower_inv(float *Lo, float *L_inv, int DIM)
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
linv_label1: for (int j = 0; j < DIM; j++)
{

	if (i < j)
		L_inv[i*DIM + j] = 0;
	else if (i == j)
		L_inv[i*DIM + j] = 1 / Lo[i*DIM + j];
	else
	{
		float sum = 0.0f;
	linv_label2: for (int k = j; k < i; k++)



		//#pragma HLS LOOP_TRIPCOUNT min=1 max=3

		//#pragma HLS PIPELINE
		sum = sum + Lo[i*DIM + k] * L_inv[k*DIM + j];

				 L_inv[i*DIM + j] = -sum;
	}
}
}
}

void Upper_inv(float *U, float *U_inv, int DIM)
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
	//#pragma HLS PIPELINE
uinv_label11: for (int j = 0; j < DIM; j++)
{
	U_inv[i*DIM + j] = 0;
}
}

		  uinv_label20: for (int i = 0; i < DIM; i++)
		  {
			  //#pragma HLS PIPELINE
			  U_inv[i*DIM + i] = (float)1 / U[i*DIM + i];
		  }

					univ_label30: for (int i = 0; i < DIM; i++)
					{
					univ_label31: for (int j = 0; j < i; j++)
					{
						//#pragma HLS PIPELINE
						//#pragma HLS LOOP_TRIPCOUNT min=1 max=3
						float sum = 0.0f;
					univ_label32: for (int k = j; k < i; k++)
					{
						//#pragma HLS PIPELINE
						//#pragma HLS LOOP_TRIPCOUNT min=1 max=3
						sum += (U[k*DIM + i] * U_inv[j*DIM + k]);
					}
								  U_inv[j*DIM + i] = -sum / U[i*DIM + i];
					}
					}
}

void matrix_mult(float *U_inv, float *L_inv, float *A_inv, int DIM)
{
	//#pragma HLS ARRAY_PARTITION variable=L_inv block factor=3 dim=1
	//#pragma HLS INLINE
	float sumTemp = 0;
MM_L1: for (int i = 0; i < DIM; i++)
{
MM_L2: for (int j = 0; j < DIM; j++)
{
	//#pragma HLS PIPELINE
	float sumFinal = 0;
MM_L3: for (int k = 0; k < DIM; k++)
{
	sumTemp = U_inv[i*DIM + k] * L_inv[k*DIM + j];
	sumFinal += sumTemp;
}
	   A_inv[i*DIM + j] = sumFinal;
}
}
}

void final_perm(float *UL_inv, float *Pe, float *A_inv, int DIM)
{
	//#pragma HLS INLINE
		/*Multiplication with P is just the row permuation*/
L1: for (int i = 0; i < DIM; i++)
{
	//#pragma HLS PIPELINE
L2: for (int j = 0; j < DIM; j++)
{
	A_inv[j*DIM + i] = UL_inv[j*DIM + (int)Pe[i]];
}
}
}

bool inverse(float *A, float *A_inv, int DIM)
{
	float *Ac = new float[DIM*DIM];
	for (int i = 0; i < DIM*DIM; i++)
	{
		Ac[i] = A[i];
	}
	//cout << "\n Entered top";
	float *Lo = new float[DIM*DIM];
	float *U = new float[DIM*DIM];
	float *L_inv = new float[DIM*DIM];
	float *U_inv = new float[DIM*DIM];
	float *UL_inv = new float[DIM*DIM];
	float *Pe = new float[DIM + 1];

	// calculate the P,A,L,U such that PA=LU
	bool singular = lup(Ac, Lo, U, Pe, DIM);

	if (singular)
	{
		cout << "SINGULAR\n";
		for (int i = 0; i < DIM; i++)
		{
			//#pragma HLS PIPELINE
			for (int j = 0; j < DIM; j++)
				A_inv[i*DIM + j] = 0;
		}
		return singular;
	}
	else
	{
		//cout << "///////////////////////////n";
		Lower_inv(Lo, L_inv, DIM);
		Upper_inv(U, U_inv, DIM);
		matrix_mult(U_inv, L_inv, UL_inv, DIM);
		final_perm(UL_inv, Pe, A_inv, DIM);


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

template<class T>
void display(T *A, int p, int l)
{
	for (int i = 0; i < p; i++)
	{
		for (int j = 0; j < l; j++)
			cout << A[i * l + j] << " ";
		cout << endl;
	}
}

template<class T>
void display_row(T *A,T*A_i, int p, int l, int j)
{
	for (int i = 0; i < p; i++)
	{
		cout << A[i * l + j] << " +i"<< A_i[i * l + j]<<"\t";
		cout << endl;
	}
}
void complexMatrixMult(float *matrix_1, float *matrix_1_i, float *matrix_2, float *matrix_2_i, int p, int l, int m, int n, float *matrix_product, float *matrix_product_i) {
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

}

void Add(float *mat1, float *mat2, int p, int l, float *mat12)
{
	for (int i = 0; i < p*l; i++)
	{
		mat12[i] = mat2[i] + mat1[i];
		//cout<<mat12[i]<<"  ";
	}
}
void neg(float *neg_inv, int m,int n, float *inv)
{
	for (int i = 0; i < m*n; i++)
	{
		inv[i] = -1 * neg_inv[i];
	}
}
void complex_inverse(float *A, float *Ai, int n, float *inv, float *inv_i)
{
	float *neg_inv_i = new float[n*n];
	float *A_i = new float[n*n];
	float *BA_i_B = new float[n*n];
	float *A_BA_i_B = new float[n*n];
	float *BA_i = new float[n*n];

	inverse(A, A_i, n);
		//display(A_i, n, n);
	//display(A,n,n);
	MatrixMult(A_i, Ai, n, n, n, n, BA_i);

	MatrixMult(Ai, BA_i, n, n, n, n, BA_i_B);
	Add(A, BA_i_B, n, n, A_BA_i_B);
	inverse(A_BA_i_B, inv, n);
	//display(A_BA_i_B,n,n); 
	//display(inv, n, n);
	MatrixMult(BA_i, inv, n, n, n, n, neg_inv_i);
	neg(neg_inv_i,n, n, inv_i);
	//display(inv_i, n, n);
}

void pseudo(float *matrix,float *matrix_i, int p, int k, float *pinv, float *pinv_i)
{


	float *matrix_mult = new float[k*k]; 	float *matrix_mult_i = new float[k*k];

	for (int i = 0; i < k; i++)
	{
		for (int j = 0; j < k; j++)
		{
			matrix_mult[i * k + j] = 0;			matrix_mult_i[i * k + j] = 0;

			for (int z = 0; z < p; z++)
			{
				matrix_mult[i * k + j] += (matrix[z * K + i] * matrix[z * K + j])-(-matrix_i[z * K + i] * matrix_i[z * K + j]);
				matrix_mult_i[i * k + j] += (-matrix_i[z * K + i] * matrix[z * K + j]) + (matrix[z * K + i] * matrix_i[z * K + j]);

			}
			//cout << matrix_mult[i * K + j] << "\t" << matrix_mult_i[i * K + j]<<"---\n";
		}
	}
	float *inv = new float[k*k]; // To store inverse 
	float *inv_i = new float[k*k]; // To store inverse 
	//cout << "\n As'*As\n";
	//display(matrix_mult, k, k); display(matrix_mult_i, k, k);
	
	
	complex_inverse(matrix_mult,matrix_mult_i,k,inv,inv_i);
	//cout << "\nINVERSE of square mat"<<k<<""<< ":\n";
	//display(inv, k, k); display(inv_i,k,k);

		
	//MatrixMult(inv, t_matrix, k, k, k, p, pinv);
	for (int i = 0; i < k; i++)
	{
		for (int j = 0; j < p; j++)
		{
			pinv[i * p + j] = 0;	pinv_i[i * p + j] = 0;
			for (int z = 0; z < k; z++)
			{
				pinv[i * p + j] += (inv[i * k + z] * matrix[j * p + z])+(inv_i[i * k + z] * matrix_i[j * p + z]);
				pinv_i[i * p + j] += ((inv_i[i * k + z] * matrix[j * p + z]) - (inv[i * k + z] * matrix_i[j * p + z]));
			}
		}
	}
	//cout << "\nThe Monroe-penrose inverse is :\n";
	//display(pinv, k, p);	display(pinv_i, k, p);



}


int max(float z[N])
{
	int maxind;
	float max = 0;
	for (int i = 0; i < N; i++)
	{
		if (z[i] > max)
		{
			max = z[i];
			maxind = i;
		}
	}
	return maxind;
}

//initializing A1 (normalised A matrix)
void INIT_A1(float A[P*N], float A_i[P*N],float A1[P*N], float A1_i[P*N])
{
	//absolute matrix of A i.e. A+Ai
	float abs_A[P*N];
	for (int i = 0; i < P*N; i++)
	{
		abs_A[i] = sqrt(A[i] * A[i] + A_i[i] * A_i[i]);
	}
	//display(abs_A,P,N);
	for (int j = 0; j < N; j++)
	{
		float b = norm(abs_A, j);
		//printf("\nNorm(A,j)= %f\t %d\n", b, j);
		for (int i = 0; i < P; i++)
		{
			A1[i*N + j] = A[i*N + j] / b;
			A1_i[i*N + j] = A_i[i*N + j] / b;
		}
	}
}

void update_Z(float A1[P*N], float A1_i[P*N],float R[P*L], float R_i[P*L],float Z[N])
{
	float product,product_i, sum,sum_i;
	// Loop for updating Z
	for (int j = 0; j < N; j++)
	{
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
		cout << endl;
		Z[j] = sqrt(sum);
	}
}

// initializing new rows in As at every iteration of k
void update_As(float A1[P*N],float As[P*K],int k, int maxind)
{
#pragma HLS PIPELINE
	for (int i = 0; i < P; i++)
	{
		As[i*K + k] = A1[i*N + maxind];

	}
}

//get Ax=As*X
void get_Ax(float As[P*N], float As_i[P*N],float X[K*L], float X_i[K*L],int k, float Ax[P*L], float Ax_i[P*L])
{
	//MatrixMult(As, X, P, k + 1, k + 1, L, Ax);
	for (int i = 0; i < P; i++)
	{
		for (int j = 0; j < L; j++)
		{
			Ax[i * L + j] = 0; Ax_i[i * L + j] = 0;
			for (int z = 0; z < (k + 1); z++)
			{
				Ax[i * L + j] += (As[i * K + z] * X[z * L + j]) - (As_i[i * K + z] * X_i[z * L + j]);
				Ax_i[i * L + j] += (As_i[i * K + z] * X[z * L + j]) + (As[i * K + z] * X_i[z * L + j]);
			}
		}
	}

}

//UPDATING RESDUE
void update_R(float R[P*L], float R_i[P*L],float Ax[P*L],float Ax_i[P*L])
{
	//UPDATING RESIDUE MATRIX
	for (int i = 0; i < P; i++)
	{
		for (int j = 0; j < L; j++)
		{
			R[i*L + j] = R[i*L + j] - Ax[i*L + j];
			R_i[i*L + j] = R_i[i*L + j] - Ax_i[i*L + j];
		}
	}
}


void OMP_sw(float A[P*N], float A_i[P*N], float R[P*L], float R_i[P*L],int out[N],int iter)
{
	float A1 [P*N];float A1_i[P*N];
	//float *A = new float[P*N];
	float As [P*K];	float As_i[P*K];

	//initializing A1 (normalised A matrix)
	INIT_A1(A,A_i,A1,A1_i);
	//display(A1, P, N);
	//display(A1_i, P, N);

	
	float Ax [P*L]; float Ax_i[P*L];

	int supp [K];
	float Z[N];
	
	int maxind;
	for (int k = 0; k < iter; k++)
	{
		float pseudoinverse[K*P];		float pseudoinverse_i[K*P];
		float X[K*L];		float X_i[K*L];


		update_Z(A1, A1_i, R, R_i, Z);
		//cout << "Z-->";
		//display(Z, N, 1);

		maxind = max(Z);
		//cout << "max--> " << maxind << endl;
		supp[k] = maxind;
		out[maxind] = 1;
		update_As(A1, As, k, maxind);
		update_As(A1_i, As_i, k, maxind);

		//display_row(As, As_i, P, K, k);

		// FINDS PSEUDO INVERSE OF THE MATRIX
		//display(As, P, k + 1);
		pseudo(As, As_i, P, k + 1, pseudoinverse, pseudoinverse_i);
		//completes the X=As\R operation
		complexMatrixMult(pseudoinverse, pseudoinverse_i, R, R_i, K, P, P, L, X, X_i);
		//cout << "As\R-->\n"; display(X, k + 1, L); display(X_i, k + 1, L);


		get_Ax(As, As_i, X, X_i, k, Ax, Ax_i);
		//cout << "Ax-->\n"; display(Ax, P, L); display(Ax_i, P, L);

		update_R(R, R_i, Ax, Ax_i);
		
		//cout << "R->\n";
		//display(R, P, L); display(R_i, P, L);
	}
	cout << "\n-----\n";
	display(supp, 1, iter);
	cout << "\nOUTPUT Bands\n";
	display(out, 1, N);
}

// Driver program
int main()
{
	/*float matrix[P][N] = { {6, -2, 2, 7,0,8},
						  {1, 0, 5, 3,1,9},
				  {-3, 1, 5, 0,2,1},
	};*/
	//float *A1 = new float[P*N];
	float *A = new float[P*N];
	float *A_i = new float[P*N];
	//float *As = new float[P*K];
	/*for (int i = 0; i < P*N; i++)
	{
		A[i] = rand() % 10 ;
		A_i[i] = rand() % 10 ;
	}*/
	float *R = new float[P*L];
	float *R_i = new float[P*L];

	/*for (int i = 0; i < P*L; i++)
	{
		R[i] = rand() % 10;
		R_i[i] = rand() % 10;
	}*/
	//float R[P*L] = { 5, -2, 2, 7,1, 0, 5, 3,-3, 1, 5, 0 };

	int *golden_out = new int[14]; int iter;
	int out[1 * N]; float count=0;
	

	FILE *fin, *fin_i, *fin2, *fin2_i, *fin3, *fin4;
	float inputsample, inputsample_i, inputsample2, inputsample2_i;
	int inputsample3, inputsample4;
#pragma warning(suppress : 4996)
	fin = fopen("data_A.dat", "r"); fin_i = fopen("data_A_i.dat", "r"); fin2 = fopen("data_R.dat", "r"); fin2_i = fopen("data_R_i.dat", "r"); fin3 = fopen("Truth.dat", "r"); fin4 = fopen("K.dat", "r");

	for (int x = 0; x < 1; x++) {
		for (int i = 0; i < N; i++)
		{
			out[i] = 0;
		}

		for (int i = 0; i < 8 * 14; i++)
		{
#pragma warning(suppress : 4996)
			fscanf(fin, "%f", &inputsample); fscanf(fin_i, "%f", &inputsample_i);
			A[i] = inputsample;		A_i[i] = inputsample_i;

		}
		for (int i = 0; i < 8 * 299; i++)
		{
#pragma warning(suppress : 4996)
			fscanf(fin2, "%f", &inputsample2); fscanf(fin2_i, "%f", &inputsample2_i);
			R[i] = inputsample2;		R_i[i] = inputsample2_i;

		}
		for (int i = 0; i < 14; i++)
		{
#pragma warning(suppress : 4996)
			fscanf(fin3, "%d", &inputsample3); 
			golden_out[i] = inputsample3;
		}
#pragma warning(suppress : 4996)
		fscanf(fin4, "%d", &inputsample4);
		iter = inputsample4;
		//display(R, P, L);
		//display(R_i, P, L);
		//converting 2D A matrix --> 1D matrix
		//Trans_2D_1D(matrix, A);
		//printf("\n");
		//display(A, P, N);
		//display(A_i, P, N);
		OMP_sw(A, A_i, R, R_i,out,iter);
		display(golden_out, 1, 14);
		int check=0;
		for (int i = 0; i < N; i++)
		{
			if (out[i] == golden_out[i])
				check++;
		}
		if (check == N)
			count++;
		else
			cout << "MISMATCH AT" << x << "th input \n";
	}
	cout << "\n\nAccuracy is " << count / 1400 * 100;
	return 0;
}
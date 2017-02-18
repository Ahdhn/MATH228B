#include "BandMatrixSolver.h"
#include <algorithm>
#include <math.h>

BandMatSolver::BandMatSolver(size_t N, size_t M1, size_t M2)
{
	my_N = N;
	my_M1 = M1;
	my_M2 = M2;

	indx = new size_t[N + 1];
	al = new double*[N + 1];
	for (size_t i = 1; i <= N; i++){
		al[i] = new double[M1 + 1];
	}
	
}

BandMatSolver::~BandMatSolver()
{
	for (size_t i = 1; i <= my_N; i++){
		free(al[i]);
	}
	free(al);
	free(indx);
}
void BandMatSolver::Solver(double**A, double*&B)
{	
	banbks(A, B);
}
void BandMatSolver::arrange(double** A)
{
	//called just once 
	size_t i, j, l;
	int mm = my_M1 + my_M2 + 1;;
	
	l = my_M1;
	//rearange the storage 	
	for (i = 1; i <= my_M1; i++){
		size_t sb = 1;
		for (j = my_M1 + 2 - i; j <= mm; j++){
			//A[i][j - 1] = A[i][j];
			A[i][sb] = A[i][j];
			sb++;
		}
		l--;
		for (j = mm - l; j <= mm; j++){
			A[i][j] = 0.0;
		}
	}

	bandec(A);
}
void BandMatSolver::bandec(double **A)
{
	size_t i, j, k, l;
	int mm = my_M1 + my_M2 + 1;;
	double dum;
	
	/*
	l = my_M1;
	//rearange the storage 	
	for (i = 1; i <= my_M1; i++){
		size_t sb = 1;
		for (j = my_M1 + 2 - i; j <= mm; j++){
			//A[i][j - 1] = A[i][j];
			A[i][sb] = A[i][j];
			sb++;
		}
		l--;
		for (j = mm - l; j <= mm; j++){ 
			A[i][j] = 0.0; 
		}
	}*/
	
	l = my_M1;
	for (k = 1; k <= my_N; k++){
		//for each row 
		dum = A[k][1];
		i = k;
		if (l < my_N){ l++; }//find the pivot element
		for (j = k + 1; j <= l; j++){ 
			if (fabs(A[j][1])>fabs(dum)){
				dum = A[j][1];
				i = j;
			}
		}
		indx[k] = i;
		if (dum == 0.0){ A[k][1] = 10E-20; }
		//matrix is algorithmically singular, but we proceed anyway with small pivot 
		if (i != k){			
			for (j = 1; j <= mm; j++){
				//interchange rows
				std::swap(A[k][j], A[i][j]);
			}
		}
		for (i = k + 1; i <= l; i++){
			//Do the elimination 
			dum = A[i][1] / A[k][1];
			al[k][i - k] = dum;
			for (j = 2; j <= mm;j++){
				A[i][j - 1] = A[i][j] - dum*A[k][j];
			}
			A[i][mm] = 0.0;
		}

	}
}
void BandMatSolver::banbks(double**A, double*&B)
{
	size_t i, k, l;	
	double dum;
	int mm = my_M1 + my_M2 + 1;
	l = my_M1;
	for (k = 1; k <= my_N; k++){
		//forward substitution, unscrambling the permuted rows as we go
		i = indx[k];
		if (i != k){ std::swap(B[k], B[i]); }
		if (l < my_N){ l++; }
		for (i = k + 1; i <= l; i++){
			B[i] -= al[k][i - k] * B[k];
		}
	}

	l = 1;
	for (i = my_N; i >= 1; i--){
		//backsubstitution
		dum = B[i];
		for (k = 2; k <= l; k++){ 
			dum -= A[i][k] * B[k + i - 1]; 
		}
		B[i] = dum / A[i][1];
		if (l < mm){ l++; }
	}	
}
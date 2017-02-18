//Copied from "Numerical Recipes in C - The Art of Scientific Computing" - 2nd Edition
//Ahmed H. Mahmoud 
//Solve the system A.X=B, where A is banded matrix of size N x N
//with M1 subdiagonal and M2 superdiagonal 
//the A is supposed to be stored in a diagonal format such that A[1..N][1...M1+M2+1]
//the diagonal are in A[1..N][M1+1], subdiagonal elements are in A[j..N][1..M1](with j>1 appropriate to
//the number of elements on each subdiagonal), superdiagonal elements are in A[1..j][M1+2...M1+M2+1] with j<N
//appropriate to the number of elements on each superdiagonal 
//and X[1..n] and B[1..n]
//The returned is solution overwrites the RHS (array B)


#ifndef _BAND_MAT_SOLVER_
#define	_BAND_MAT_SOLVER_

class BandMatSolver
{
public:
	BandMatSolver(size_t, size_t, size_t);
	~BandMatSolver();

	//Call this once to re-arrange the matrix A for subsequenct solve
	//if matrix A is used more than once, then you just need to call this once and then call 
	//Solver() with different B's
	void arrange(double** A); 

	//Call this to solve	
	void Solver(double**A, double*&B);
	


private:
	size_t*indx;
	double**al;
	size_t my_N, my_M1, my_M2;
	
	void bandec(double **A);
	void banbks(double**A, double*&B);
};



#endif 
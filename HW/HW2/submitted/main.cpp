//Solving problem 2 in homework two - MATH 228B - Winter 2017
//By Ahmed H. Mahmoud
//The code solves the problems, and produces raw data file for matlab to plot 

//#include "Plotter.h"
#include "BandMatrixSolver.h"

#include <iostream>
#include <fstream>
#include <string.h> 
#define _USE_MATH_DEFINES
#include<math.h>
#include <algorithm>
#include <time.h>

#define _tol 10E-20


using namespace std;


size_t _SCHEME;//scheme's flag



inline void PrintOut(string outputfilename,double***U,size_t n_x,size_t n_y,double del_x,double del_y,size_t time)
{
	//print out the U at certain time
	ofstream file(outputfilename,ios::out);
	file.precision(10);
	
	file << n_x << " " << n_y << endl;

	for (size_t j = 0; j < n_y; j++){ 
		double y = j*del_y;
		for (size_t i = 0; i < n_x; i++){ 
			double x = i*del_x;
			file << U[time][i][j] <<"    ";
		}
		file << endl;
	}
}
inline void PrintOutTransAtPosition(double x,double y,double***U,size_t n_x,size_t n_y,double del_x,double del_y,double del_t,size_t steady_time)
{
	//print out the transient solution (various times) for a certain position (x,y)
	size_t i, j, t;
	i = size_t(x / del_x);
	j = size_t(y / del_y);
	fstream file("Trans_Solution_fixed_point.txt", ios::out);
	file.precision(5);
	file << "(" << x << " " << y << ")" << endl;
	file << "Time U" << endl;

	for (t = 0; t<steady_time; t++){
		file << t*del_t << " " << U[t][i][j] << endl;
	}
	file.close();

}

inline void ImposeBoundaryCondition(size_t n_x,double del_x, size_t n_y,double del_y, double***U, size_t time_size)
{
	//set the initial and boundary conditions 
	
	//at time>=0, U=0 @ 0 <= x <=X && y=0
	//at time>=0, U=0 @ 0 <  y < Y && x=0
	//at time>=0, U=0 @ 0 <= x <=X && y=Y
	//at time>=0, U=0 @ 0 <  y < Y && x=X

	//at time =0, T=exp(-100(x-0.3)^2 + (y-0.4)^2)) @ 0< x< X && 0<y<Y
	

	//1) initial conditions 
	for (size_t i = 1; i < n_x - 1; i++){ 
		for (size_t j = 1; j < n_y - 1; j++){ 
			double x = double(i)*del_x;
			double y = double(j)*del_y;
			double xx = (x - 0.3)*(x - 0.3);
			double yy = (y - 0.4)*(y - 0.4);
			U[0][i][j] = exp(-100 * (xx + yy));
		}
	}

	//2) boundary conditions 
	//at time>=0, U=0 at 0 <= x <=X && y=0
	for (size_t t = 0; t < time_size; t++){
		for (size_t i = 0; i < n_x; i++){ 
			U[t][i][0]=0;
		}
	}

	//at time>=0, U=0 at 0 <  y < Y && x=0
	for (size_t t = 0; t < time_size; t++){
		for (size_t j = 1; j < n_y - 1; j++){ 
			U[t][0][j] = 0;
		}
	}

	//at time>=0, U=0 at 0 <= x <=X && y=Y
	for (size_t t = 0; t < time_size; t++){
		for (size_t i = 0; i < n_x; i++){ 
			U[t][i][n_y - 1] = 0;
		}
	}	

	//at time>=0, U=0 at 0 <  y < Y && x=X
	for (size_t t = 0; t < time_size; t++){
		for (size_t j = 1; j < n_y - 1; j++){ 
			U[t][n_x - 1][j] = 0;
		}
	}
}

inline void ForwardEuler(double***U, size_t n_x, size_t n_y, double alpha, double del_x, double del_y, double del_t, size_t max_t, double&my_time_mil_sec)
{
	//solving using forward Euler explicit Scheme
	my_time_mil_sec = 0;
	double cof_x, cof_y;//coefficients
	cof_x = (alpha*del_t) / (del_x*del_x);
	cof_y = (alpha*del_t) / (del_y*del_y);

	if (cof_x + cof_y > 0.5 - _tol){
		double tona;
		cout << "\n WARNING at ForwardEuler()!! The solution produced won't be stable since" << endl;
		cout << "(alpha*delta_t)/(delta_x^2) + (alpha*delta_t)/(delta_y^2)= " << cof_x + cof_y << ", thus it's > 0.5" << endl;
		cout << "\nPress any ket to exit" << endl;
		cin >> tona;
		exit(0); 
	}
		
	size_t my_current_time;
	size_t my_prev_time;

	 
	clock_t t1, t2;
	t1 = clock();//get time 

	for (size_t t = 1; t < max_t; t++){ 
		if (t % 2 == 1){
			my_current_time = 1;
			my_prev_time = 0;
		}
		else{
			my_current_time = 0;
			my_prev_time = 1;
		}
		
		for (size_t i = 1; i < n_x - 1; i++){ 
			for (size_t j = 1; j < n_y - 1; j++){ 
				U[my_current_time][i][j] = U[my_prev_time][i][j] +
					                       cof_x*(U[my_prev_time][i + 1][j] - 2 * U[my_prev_time][i][j] + U[my_prev_time][i - 1][j]) +
					                       cof_y*(U[my_prev_time][i][j + 1] - 2 * U[my_prev_time][i][j] + U[my_prev_time][i][j - 1]);
			}
		}	
	}
	t2 = clock();//get time 

	my_time_mil_sec = (double(t2) - double(t1)) / double(CLOCKS_PER_SEC * 1000.0);

}

inline void CrankNicolson(double***U, size_t n_x, size_t n_y, double del_x, double del_y, double del_t, size_t max_t, double&my_time_mil_sec)
{
	//solve for U using CN implicit scheme  
	//at each time step we get a penta-diagonal system to be solved
	my_time_mil_sec = 0;

	if (n_x != n_y){
		size_t tona;
		cout << "\n WARNING at CrankNicolson()!! Can not work when spatial discretization steps in X is different than Y " << endl;		
		cout << "\nPress any ket to exit" << endl;
		cin >> tona;
		exit(0);
	}
	
	double r(del_t / del_x*del_x);
	double alpha(1 + 2.0*r), beta(-r / 2.0), gamma(1 - 2.0 * r);
	
	//matrix passed to the solver (full of constants)
	//length = number of unknows = (n_x-2)*(n_y-2) [numbering starts from 1]
	//width is the bandwidth of the banded matrix (idealy it's 5, but we also pass some zeros)
	//width = (n_x-2)+(n_y-2) + 1	

	size_t M1(n_x - 2), M2(n_y - 2), N((n_x - 2)*(n_y - 2));
	double**penta_matrix = new double*[N+1];
	for (size_t i = 1; i <= N; i++){ 
		penta_matrix[i] = new double[M1 + M2 + 2];
		//we initiate the solution with zeros so we don't worry about
		//fill in the zeros when they are actually needed 
		for (size_t j = 1; j <= M1 + M2 + 1; j++){ penta_matrix[i][j] = 0; }
	}
	double*B = new double[N+1];//the LHS input vector (solution overwrites this vector)
	
	//**********
	//fill in the banded matrix 
	for (size_t i = 1; i <= N; i++){
		//the diagonal
		penta_matrix[i][M1 + 1] = alpha;
		if (i != 1){
			//the 1st subdiagonal
			if ((i - 1) % M1 == 0){
				penta_matrix[i][M1] = 0;
			}
			else{
				penta_matrix[i][M1] = beta;
			}
		}
		if (i != N){
			//the 1st superdiagonal
			if (i % M1 == 0){ 
				penta_matrix[i][M1 + 2] = 0;
			}
			else{
				penta_matrix[i][M1 + 2] = beta;
			}
		}
	}
	
	for (size_t i = M1 + 1; i <= N; i++){ 
		//the outer subdiagonal		
		penta_matrix[i][1] = beta;		
	}
	for (size_t i = 1; i <= N - M2 - 1; i++){
		//the outer superdiagonal
		penta_matrix[i][M1 + M2 + 1] = beta;
	}
	//**********


	BandMatSolver pentaMatSol(N, M1, M2);
	pentaMatSol.arrange(penta_matrix);

	
	size_t my_current_time = 1;
	size_t my_prev_time = 0;

	for (size_t t = 1; t < max_t; t++){
		clock_t t1, t2;
		t1 = clock();

		if (t % 2 == 1){
			my_current_time = 1;
			my_prev_time = 0;
		}
		else{
			my_current_time = 0;
			my_prev_time = 1;
		}
		
		//fill in the RHS vector B	
		for (size_t i = 1; i <= M1; i++){
			for (size_t j = 1; j <= M2; j++){
				size_t id = (i - 1)*(M1)+j;
				B[id] = gamma*U[my_prev_time][i][j] - beta*(U[my_prev_time][i + 1][j] + U[my_prev_time][i - 1][j] + U[my_prev_time][i][j + 1] + U[my_prev_time][i][j - 1]);

				if (i == 1){//i-1, j
					B[id] -= beta * U[my_current_time][0][j];
				}
				if (j == 1){//i, j-1
					B[id] -= beta * U[my_current_time][i][0];
				}
				if (i == M1){//i+1, j
					B[id] -= beta * U[my_current_time][n_x - 1][j];
				}
				if (j == M2){//i, j+1
					B[id] -= beta * U[my_current_time][i][n_y - 1];
				}

			}
		}

		//solve the matrix for this time step
		pentaMatSol.Solver(penta_matrix, B);
		t2 = clock();

		//only record time for computation (not for data movement)
		my_time_mil_sec += (double(t2) - double(t1)) / double(CLOCKS_PER_SEC * 1000.0);

		//copy the solution to main U
		for (size_t i = 1; i <= M1; i++){
			for (size_t j = 1; j <= M2; j++){ 
				size_t id = (i - 1)*(M1)+j;
				U[my_current_time][i][j] = B[id];
			}
		}
	}

	
	for (size_t i = 1; i <= N; i++){
		free(penta_matrix[i]);
	}
	delete[]penta_matrix;
	delete[] B;
		
}
void TestPentSolver()
{
	/*
	Matrix A
	1 8 0 6 0 0 0 0 0;
	1 8 4 0 4 0 0 0 0;
	0 8 1 0 0 1 0 0 0;
	1 0 0 6 4 0 5 0 0;
	0 8 0 6 4 5 0 5 0;
	0 0 4 0 3 5 0 0 1;
	0 0 0 6 0 0 4 8 0;
	0 0 0 0 1 0 4 8 9;
	0 0 0 0 0 4 0 9 4;

	Vector B
	30
	1
	4
	6
	96
	6
	31
	88
	20
	*/

	size_t M1(5), M2(5), N(25);
	

	double**mat = new double*[N + 1];
	for (size_t i = 1; i <= N; i++){ 
		mat[i] = new double[M1 + M2 + 2];
	}
	double*B = new double[N + 1];

	mat[1][1] = 0;	                        mat[1][2] = 0 ;	mat[1][3] = 0 ;	mat[1][4] = 0 ;	mat[1][5] = 0;	                        mat[1][6] = 1.0002;	    mat[1][7] = -5.0000000000000002e-005;	mat[1][8] = 0;	mat[1][9] = 0;	mat[1][10] = 0;	mat[1][11] = -5.0000000000000002e-005;
	mat[2][1] = 0;	                        mat[2][2] = 0 ;	mat[2][3] = 0 ;	mat[2][4] = 0 ;	mat[2][5] = -5.0000000000000002e-005;	mat[2][6] = 1.0002;	    mat[2][7] = -5.0000000000000002e-005;	mat[2][8] = 0;	mat[2][9] = 0;	mat[2][10] = 0;	mat[2][11] = -5.0000000000000002e-005;
	mat[3][1] = 0;	                        mat[3][2] = 0 ;	mat[3][3] = 0 ;	mat[3][4] = 0 ;	mat[3][5] = -5.0000000000000002e-005;	mat[3][6] = 1.0002;	    mat[3][7] = -5.0000000000000002e-005;	mat[3][8] = 0;	mat[3][9] = 0;	mat[3][10] = 0;	mat[3][11] = -5.0000000000000002e-005;
	mat[4][1] = 0;	                        mat[4][2] = 0 ;	mat[4][3] = 0 ;	mat[4][4] = 0 ;	mat[4][5] = -5.0000000000000002e-005;	mat[4][6] = 1.0002;	    mat[4][7] = -5.0000000000000002e-005;	mat[4][8] = 0;	mat[4][9] = 0;	mat[4][10] = 0;	mat[4][11] = -5.0000000000000002e-005;
	mat[5][1] = 0;	                        mat[5][2] = 0 ;	mat[5][3] = 0 ;	mat[5][4] = 0 ;	mat[5][5] = -5.0000000000000002e-005;	mat[5][6] = 1.0002;	    mat[5][7] = 0;	mat[5][8] = 0;	mat[5][9] = 0;	mat[5][10] = 0;	mat[5][11] = -5.0000000000000002e-005;
	mat[6][1] = -5.0000000000000002e-005;	mat[6][2] = 0 ;	mat[6][3] = 0 ;	mat[6][4] = 0 ;	mat[6][5] = 0;	                        mat[6][6] = 1.0002;	    mat[6][7] = -5.0000000000000002e-005;	mat[6][8] = 0;	mat[6][9] = 0;	mat[6][10] = 0;	mat[6][11] = -5.0000000000000002e-005;
	mat[7][1] = -5.0000000000000002e-005;	mat[7][2] = 0 ;	mat[7][3] = 0 ;	mat[7][4] = 0 ;	mat[7][5] = -5.0000000000000002e-005;	mat[7][6] = 1.0002;	    mat[7][7] = -5.0000000000000002e-005;	mat[7][8] = 0;	mat[7][9] = 0;	mat[7][10] = 0;	mat[7][11] = -5.0000000000000002e-005;
	mat[8][1] = -5.0000000000000002e-005;	mat[8][2] = 0 ;	mat[8][3] = 0 ;	mat[8][4] = 0 ;	mat[8][5] = -5.0000000000000002e-005;	mat[8][6] = 1.0002;	    mat[8][7] = -5.0000000000000002e-005;	mat[8][8] = 0;	mat[8][9] = 0;	mat[8][10] = 0;	mat[8][11] = -5.0000000000000002e-005;
	mat[9][1] = -5.0000000000000002e-005;	mat[9][2] = 0 ;	mat[9][3] = 0 ;	mat[9][4] = 0 ;	mat[9][5] = -5.0000000000000002e-005;	mat[9][6] = 1.0002;	    mat[9][7] = -5.0000000000000002e-005;	mat[9][8] = 0;	mat[9][9] = 0;	mat[9][10] = 0;	mat[9][11] = -5.0000000000000002e-005;
	mat[10][1] = -5.0000000000000002e-005;	mat[10][2] = 0;	mat[10][3] = 0;	mat[10][4] = 0;	mat[10][5] = -5.0000000000000002e-005;	mat[10][6] = 1.0002;	mat[10][7] = 0;	mat[10][8] = 0;	mat[10][9] = 0;	mat[10][10] = 0;	mat[10][11] = -5.0000000000000002e-005;
	mat[11][1] = -5.0000000000000002e-005;	mat[11][2] = 0;	mat[11][3] = 0;	mat[11][4] = 0;	mat[11][5] = 0;	                        mat[11][6] = 1.0002;	mat[11][7] = -5.0000000000000002e-005;	mat[11][8] = 0;	mat[11][9] = 0;	mat[11][10] = 0;	mat[11][11] = -5.0000000000000002e-005;
	mat[12][1] = -5.0000000000000002e-005;	mat[12][2] = 0;	mat[12][3] = 0;	mat[12][4] = 0;	mat[12][5] = -5.0000000000000002e-005;	mat[12][6] = 1.0002;	mat[12][7] = -5.0000000000000002e-005;	mat[12][8] = 0;	mat[12][9] = 0;	mat[12][10] = 0;	mat[12][11] = -5.0000000000000002e-005;
	mat[13][1] = -5.0000000000000002e-005;	mat[13][2] = 0;	mat[13][3] = 0;	mat[13][4] = 0;	mat[13][5] = -5.0000000000000002e-005;	mat[13][6] = 1.0002;	mat[13][7] = -5.0000000000000002e-005;	mat[13][8] = 0;	mat[13][9] = 0;	mat[13][10] = 0;	mat[13][11] = -5.0000000000000002e-005;
	mat[14][1] = -5.0000000000000002e-005;	mat[14][2] = 0;	mat[14][3] = 0;	mat[14][4] = 0;	mat[14][5] = -5.0000000000000002e-005;	mat[14][6] = 1.0002;	mat[14][7] = -5.0000000000000002e-005;	mat[14][8] = 0;	mat[14][9] = 0;	mat[14][10] = 0;	mat[14][11] = -5.0000000000000002e-005;
	mat[15][1] = -5.0000000000000002e-005;	mat[15][2] = 0;	mat[15][3] = 0;	mat[15][4] = 0;	mat[15][5] = -5.0000000000000002e-005;	mat[15][6] = 1.0002;	mat[15][7] = 0;	mat[15][8] = 0;	mat[15][9] = 0;	mat[15][10] = 0;	mat[15][11] = -5.0000000000000002e-005;
	mat[16][1] = -5.0000000000000002e-005;	mat[16][2] = 0;	mat[16][3] = 0;	mat[16][4] = 0;	mat[16][5] = 0;	                        mat[16][6] = 1.0002;	mat[16][7] = -5.0000000000000002e-005;	mat[16][8] = 0;	mat[16][9] = 0;	mat[16][10] = 0;	mat[16][11] = -5.0000000000000002e-005;
	mat[17][1] = -5.0000000000000002e-005;	mat[17][2] = 0;	mat[17][3] = 0;	mat[17][4] = 0;	mat[17][5] = -5.0000000000000002e-005;	mat[17][6] = 1.0002;	mat[17][7] = -5.0000000000000002e-005;	mat[17][8] = 0;	mat[17][9] = 0;	mat[17][10] = 0;	mat[17][11] = -5.0000000000000002e-005;
	mat[18][1] = -5.0000000000000002e-005;	mat[18][2] = 0;	mat[18][3] = 0;	mat[18][4] = 0;	mat[18][5] = -5.0000000000000002e-005;	mat[18][6] = 1.0002;	mat[18][7] = -5.0000000000000002e-005;	mat[18][8] = 0;	mat[18][9] = 0;	mat[18][10] = 0;	mat[18][11] = -5.0000000000000002e-005;
	mat[19][1] = -5.0000000000000002e-005;	mat[19][2] = 0;	mat[19][3] = 0;	mat[19][4] = 0;	mat[19][5] = -5.0000000000000002e-005;	mat[19][6] = 1.0002;	mat[19][7] = -5.0000000000000002e-005;	mat[19][8] = 0;	mat[19][9] = 0;	mat[19][10] = 0;	mat[19][11] = -5.0000000000000002e-005;
	mat[20][1] = -5.0000000000000002e-005;	mat[20][2] = 0;	mat[20][3] = 0;	mat[20][4] = 0;	mat[20][5] = -5.0000000000000002e-005;	mat[20][6] = 1.0002;	mat[20][7] = 0;	mat[20][8] = 0;	mat[20][9] = 0;	mat[20][10] = 0;	mat[20][11] = 0;
	mat[21][1] = -5.0000000000000002e-005;	mat[21][2] = 0;	mat[21][3] = 0;	mat[21][4] = 0;	mat[21][5] = 0;	                        mat[21][6] = 1.0002;	mat[21][7] = -5.0000000000000002e-005;	mat[21][8] = 0;	mat[21][9] = 0;	mat[21][10] = 0;	mat[21][11] = 0;
	mat[22][1] = -5.0000000000000002e-005;	mat[22][2] = 0;	mat[22][3] = 0;	mat[22][4] = 0;	mat[22][5] = -5.0000000000000002e-005;	mat[22][6] = 1.0002;	mat[22][7] = -5.0000000000000002e-005;	mat[22][8] = 0;	mat[22][9] = 0;	mat[22][10] = 0;	mat[22][11] = 0;
	mat[23][1] = -5.0000000000000002e-005;	mat[23][2] = 0;	mat[23][3] = 0;	mat[23][4] = 0;	mat[23][5] = -5.0000000000000002e-005;	mat[23][6] = 1.0002;	mat[23][7] = -5.0000000000000002e-005;	mat[23][8] = 0;	mat[23][9] = 0;	mat[23][10] = 0;	mat[23][11] = 0;
	mat[24][1] = -5.0000000000000002e-005;	mat[24][2] = 0;	mat[24][3] = 0;	mat[24][4] = 0;	mat[24][5] = -5.0000000000000002e-005;	mat[24][6] = 1.0002;	mat[24][7] = -5.0000000000000002e-005;	mat[24][8] = 0;	mat[24][9] = 0;	mat[24][10] = 0;	mat[24][11] = 0;
	mat[25][1] = -5.0000000000000002e-005;	mat[25][2] = 0;	mat[25][3] = 0;	mat[25][4] = 0;	mat[25][5] = -5.0000000000000002e-005;	mat[25][6] = 1.0002;	mat[25][7] = 0;	mat[25][8] = 0;	mat[25][9] = 0;	mat[25][10] = 0;	mat[25][11] = 0;


	B[1] = 0.000735644;
	B[2] = 0.108378183;
	B[3] = 0.062185974;
	B[4] = 0.000141031;
	B[5] = 8.08E-09;
	B[6] = 0.003893875;
	B[7] = 0.573661329;
	B[8] = 0.329159319;
	B[9] = 0.000746499;
	B[10] = 4.28E-08;
	B[11] = 7.99E-05;
	B[12] = 0.011770308;
	B[13] = 0.006753647;
	B[14] = 1.53E-05;
	B[15] = 8.76E-10;
	B[16] = 1.03E-08;
	B[17] = 1.52E-06;
	B[18] = 8.70E-07;
	B[19] = 1.96E-09;
	B[20] = 7.57E-14;
	B[21] = 3.15E-13;
	B[22] = 4.67E-11;
	B[23] = 2.68E-11;
	B[24] = 5.95E-14;
	B[25] = 5.28E-19;	

	BandMatSolver pentaMatSol(N, M1, M2);
	pentaMatSol.arrange(mat);
	pentaMatSol.Solver(mat, B);
	fstream file_B("B.txt", ios::out);
	file_B.precision(30);
	for (size_t i = 1; i <= N; i++){
		file_B << B[i] << endl;
	}
	file_B.close();
}
int main (int argc, char *argv[])
{
	//TestPentSolver();

	double***U;//solution for n time to all grid points	
	double X, Y, total_time, del_t, del_x, del_y;
	size_t n_x, n_y, max_t, steady_time; 


	cout << "Solving 2 in homework two - MATH 228B - Winter 2017" << endl;
	_SCHEME=0;	
	
	cout<<"\n******** SCHEME ********"<<endl;
	while (_SCHEME != 1 && _SCHEME != 2){ 
		cout<<"\nPlease enter number 1 for forward Euler or 2 for Crank-Nicolson: "<<endl;		
		cin>>_SCHEME;
		if (_SCHEME != 1 && _SCHEME != 2){ 
			cout<<"Invalid Entry!!!"<<endl;
		}
	}
	
	//GEOMETRY & BOUNDARY CONDITIONS
		
	total_time = 1.0; //total time of simulation 
	X = 1.0; //length in x-direction
	Y = 1.0;//width in y-direction 

		
	U = new double**[2];//we only save solution at t and t-1
	
	del_x = 1;
	del_y = 1;

	if (_SCHEME){
		del_x = 0.5;
		del_y = 0.5;
	}

	//at each itteration, get new space step
	//solve and report time 

	string file_name = (_SCHEME == 1) ? "timing_FE.txt" : "timing_CN.txt";

	fstream timing_file(file_name, ios::out);
	timing_file.precision(30);
	
	timing_file << "ittr     del_t      del_x      time(MSec)" << endl;
	

	for (size_t ittr = 1; ittr <= 10; ittr++){ 

		//refine space
		del_x /= 2.0;
		del_y = del_x;

		//decide the time step
		if (_SCHEME == 2){
			//for CN, use fixed time step at all refine levels 
			del_t = 0.01;	
		}
		else{
			//for forward euler, time step is just below the stability limit 
			del_t = ((del_x*del_x) / 4.0) - 2.0*_tol;
		}
		
		n_x = size_t(ceil(X / del_x) + 1);
		n_y = size_t(ceil(Y / del_y) + 1);
		max_t = size_t(ceil(total_time / del_t) + 1);

		for (size_t t = 0; t < 2; t++){
			U[t] = new double*[n_x];
			for (size_t i = 0; i < n_x; i++){
				U[t][i] = new double[n_y];
			}
		}
		ImposeBoundaryCondition(n_x, del_x, n_y, del_y, U, 2);


		double time_mil_sec;
		if (_SCHEME == 1){
			cout << "max_t= " << max_t << " >****<  n_x=" << n_x << endl;
			ForwardEuler(U, n_x, n_y, 1.0, del_x, del_y, del_t, max_t, time_mil_sec);
			timing_file <<ittr <<"        " <<del_t << "        " << del_x << "        " << time_mil_sec << endl;
		}
		else{
			cout << "max_t= " << max_t << " >****<  n_x=" << n_x << endl;
			CrankNicolson(U, n_x, n_y, del_x, del_y, del_t, max_t, time_mil_sec);
			timing_file << ittr << "        " << del_t << "        " << del_x << "        " << time_mil_sec << endl;
		}
		for (size_t t = 0; t < 2; t++){
			for (size_t i = 0; i < n_x; i++){
				delete[] U[t][i];
			}
			delete[] U[t];
		}
	}
	


	double dummy;
	cout << "\n Press any key to quit" << endl;
	cin >> dummy;
	return 0;


}

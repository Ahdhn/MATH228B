//Solving problem 2 in homework three - MATH 228B - Winter 2017
//By Ahmed H. Mahmoud
//The code solves the problems, and produces raw data file for matlab to plot 

#include "Plotter.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h> 
#include <algorithm>
#define _USE_MATH_DEFINES
#include<math.h>
using namespace std;



Plotter _plot;//for plotting anywhere

inline void PrintOut(string outputfilename,double***U,size_t n_x,size_t n_y,double del_x,double del_y,size_t time)
{
	//print out raw data for matlab to draw 
	ofstream file(outputfilename,ios::out);
	file.precision(10);		
	file << n_x << "  " << n_y << endl;

	for (size_t j = 0; j < n_y; j++){ 		
		for (size_t i = 0; i < n_x; i++){ 			
			file<<U[time][i][j]<<"   ";			
		}
		file<<endl;
	}
}

inline void ImposeBoundaryCondition(size_t n_x, double del_x, size_t n_y, double del_y, double***U, double***W, size_t time_size, bool part_1)
{
	//set the initial	
	//if (part_1) -> at time =0, U=exp(-100(x)^2 + (y)^2)), W=0 @ 0< x< X && 0<y<Y
	//else -> U=1-2x, W=0.05

	//1) initial conditions 	
	for (size_t i = 0; i < n_x; i++){ 
		for (size_t j = 0; j < n_y; j++){ 
			double x = double(i)*del_x;
			double y = double(j)*del_y;			
			if (part_1){
				U[0][i][j] = exp(-100 * (x*x + y*y));
				W[0][i][j] = 0;
			}
			else{
				U[0][i][j] = 1.0 - 2.0*x;
				W[0][i][j] = 0.05*y;
			}
		}
	}
}

void SolveTridiagonal(double*x, size_t X, double sub, double alpha_1st_last, double alpha_mid, double super, double*helper)
{
	/*Modified from (link below) such that it handles only systems where the subdiagonal (sub) and superdiagonal (super)
	//are constants, and the main diagonal is a constant value also such that the first and last entries (alpha_1st_last) differ
	//than the values in between (alpha_mid), and all the values in between are the same 
	https://en.wikibooks.org/wiki/Algorithm_Implementation/Linear_Algebra/Tridiagonal_matrix_algorithm
	solves Ax = v where A is a tridiagonal matrix consisting of vectors a, b, c
	x - initially contains the input vector v, and returns the solution x. indexed from 0 to X - 1 inclusive
	X - number of equations (length of vector x)
	sub/a - subdiagonal (means it is the diagonal below the main diagonal), indexed from 1 to X - 1 inclusive (Modified: we take it as single scalar value)
	alpha/b - the main diagonal, indexed from 0 to X - 1 inclusive
	super/c - superdiagonal (means it is the diagonal above the main diagonal), indexed from 0 to X - 2 inclusive (Modified: we take it as single scalar value)

	Note: contents of input vector c will be modified, making this a one-time-use function (scratch space can be allocated instead for this purpose to make it reusable)
	Note 2: We don't check for diagonal dominance, etc.; this is not guaranteed stable
	*/

	helper[0] = super / alpha_1st_last;
	x[0] /= alpha_1st_last;

	/* loop from 1 to X - 1 inclusive, performing the forward sweep */
	for (size_t ix = 1; ix < X; ix++) {
		double m = 1.0 / (((ix == X - 1) ? alpha_1st_last : alpha_mid) - sub * helper[ix - 1]);
		helper[ix] = super * m;
		x[ix] = (x[ix] - sub * x[ix - 1]) * m;
	}

	/* loop from X - 2 to 0 inclusive (safely testing loop condition for an unsigned integer), to perform the back substitution */
	for (int ix = X - 2; ix >= 0; ix--){ 
		x[ix] -= helper[ix] * x[ix + 1];
	}
}

void ADI_singleStep(double***U, size_t n_x, size_t n_y, double del_x, double del_y, double del_t, size_t max_t, double*helper_container, double*b, double D)
{
	//solving using ADI for single time step
	//start with x-sweep then y-sweep
	//solving the tri-dia matrix using LU-decomposition method
	//http://www.math.buffalo.edu/~pitman/courses/mth437/na2/node2.html
	//at each sweep level, we only save two time steps
	//1) for x-sweep, we have n-time step at U[0] and save (n+0.5)-time step at U[1]
	//2) for y-sweep, we have (n+0.5)-time step at U[1] and save (n+1)-time step at U[0] 
	
	if (n_x != n_y){
		std::cout << "Error (1) ADI(). n_x != n_y" << endl;				
		system("pause");
		exit(0);	
	}

	double r,alfa, beta, gamma; 
	
	r = (D*del_t) / (2.0*del_x*del_x);
	alfa = 1.0 + r;
	beta = 1.0 + 2.0*r;
	gamma = 1.0 - 2.0*r;
		
		
	//for (size_t t = 1; t < max_t; t++){ 
	//x-sweep		
	for (size_t j = 0; j < n_y; j++){ //here we solve from 0 -->n_y-1
		for (size_t i = 0; i < n_x; i++){  //setting up the right hand side 
			b[i] = gamma*U[0][i][j];
			if (j == 0){
				b[i] += r*U[0][i][j];
			}
			else{
				b[i] += r*U[0][i][j - 1];
			}
			if (j == n_y - 1){
				b[i] += r*U[0][i][j];
			}
			else{
				b[i] += r*U[0][i][j + 1];
			}
		}
		SolveTridiagonal(b, n_x, -r, alfa, beta, -r, helper_container);
		for (size_t id = 0; id < n_x; id++){ U[1][id][j] = b[id]; }//copy solution (this is wrongly included in time computation)			
	}

	//y-sweep
	for (size_t i = 0; i < n_x; i++){ //here we solve from 0 -->n_x-1
		for (size_t j = 0; j < n_y; j++){ //setting up the right hand side 
			//here we store the RHS in place of solution where it is manipulated
			//in SolveTridiagonal
			U[0][i][j] = gamma*U[1][i][j];
			if (i == 0){
				U[0][i][j] += r*U[1][i][j];
			}
			else{
				U[0][i][j] += r*U[1][i - 1][j];
			}
			if (i == n_x - 1){
				U[0][i][j] += r*U[1][i][j];
			}
			else{
				U[0][i][j] += r*U[1][i + 1][j];
			}
		}

		SolveTridiagonal(U[0][i], n_y, -r, alfa, beta, -r, helper_container);
	}
	//}
}

inline double Reaction(double v, double w, double a, double I)
{
	//return the reaction part of FitzHugh-Nagumo equation (dv/dt)
	return (a - v)*(v - 1)*v - w + I;
}
inline double dW_dt(double v, double w, double eps, double gamma)
{
	//return dW/dt = eps(v-gamma*w)
	return eps*(v - gamma*w);
}
void RK4_reaction_W(size_t n_x, size_t n_y, double del_t, double **U, double **W, double a, double I, double eps, double gamma)
{
	//the coupled reaction and w using RK4 
	//the update happens in place 
	for (size_t i = 0; i < n_x; i++){
		for (size_t j = 0; j < n_y; j++){
			double k1v = del_t*Reaction(U[i][j], W[i][j], a, I);
			double k1w = del_t*dW_dt(U[i][j], W[i][j], eps, gamma);

			double k2v = del_t*Reaction(U[i][j] + 0.5*k1v, W[i][j] + 0.5*k1w, a, I);
			double k2w = del_t*dW_dt(U[i][j] + 0.5*k1v, W[i][j] + 0.5*k1w, eps, gamma);

			double k3v = del_t*Reaction(U[i][j] + 0.5*k2v, W[i][j] + 0.5*k2w, a, I);
			double k3w = del_t*dW_dt(U[i][j] + 0.5*k2v, W[i][j] + 0.5*k2w, eps, gamma);

			double k4v = del_t*Reaction(U[i][j] + k3v, W[i][j] + k3w, a, I); 
			double k4w = del_t*dW_dt(U[i][j] + k3v, W[i][j] + k3w, eps, gamma); 

			U[i][j] = U[i][j] + (1.0 / 6.0)*(k1v + 2.0*k2v + 2.0*k3v + k4v);
			W[i][j] = W[i][j] + (1.0 / 6.0)*(k1w + 2.0*k2w + 2.0*k3w + k4w);
		}
	}
}

int main (int argc, char *argv[])
{
	double***U,***W;//solution at time step for all grid points for v and w	
	double X, Y, del_t, del_x, del_y, total_time;
	size_t n_x, n_y, max_t, n; 
	double*helper_container, *b;//helper array for tri-diag matrix solution 	

	//user input values	
	std::cout << "Solving problem 2 in homework three - MATH 228B - Winter 2017" << endl;

	double a(0.1), gam(2.0), eps(0.005), D(5 * 10E-5),I(0.0);//equation parameters 


	total_time = 600.0;//total time 	
	X = 1.0;//length in x-direction
	Y = X;//length in y-direction
	
	del_x = 0.01;
	del_t = 0.1;

	del_y = del_x;//same in x and y

	n_x = size_t(ceil(X / del_x) + 1); //gird num points 
	n_y = n_x;

	max_t = size_t(ceil(total_time / del_t) + 1);
	n = n_x*n_y; //total spacing points 

	
	U = new double**[2];//we only need to store two time steps
	                    //either [n] and [n+(1/2)] or [n+(1/2)] and [n+1]
	W = new double**[2];
	
	//initlizing containers
	helper_container = new double[n - 1];
	b = new double[n - 1];

	for (size_t t = 0; t < 2; t++){
		U[t] = new double*[n_x];
		W[t] = new double*[n_x];
		for (size_t i = 0; i < n_x; i++){
			U[t][i] = new double[n_y];
			W[t][i] = new double[n_y];
			for (size_t j = 0; j < n_y; j++){
				U[t][i][j] = 0;
				W[t][i][j] = 0;
			}
		}
	}

	ImposeBoundaryCondition(n_x, del_x, n_y, del_y, U, W, 2, 0); 
	


	//Using Strang Splitting 

	//_plot.ColorIsoCountoring("strang.ps", del_x, n_x, del_y, n_y, U[0]);
	string fname = "p2_part3_" + to_string(0);
	PrintOut(fname + ".txt", U, n_x, n_y, del_x, del_y, 0);
	
	
	for (size_t t = 1; t < max_t; t++){		
		//1) Solve for reaction and w with time step del_t/2 (coupled system of ODE for U and w)
		RK4_reaction_W(n_x, n_y, 0.5*del_t, U[0], W[0], a, I, eps, gam);

		//2) Solve for diffusion with time step del_t
		ADI_singleStep(U, n_x, n_y, del_x, del_y, del_t, max_t, helper_container, b, D);

		//3) Solve for reaction and w with time step del_t/2 (coupled system of ODE for U and w)
		RK4_reaction_W(n_x, n_y, 0.5*del_t, U[0], W[0], a, I, eps, gam);

		
		if (t % (max_t/20) == 0){
			//_plot.ColorIsoCountoring("strang.ps", del_x, n_x, del_y, n_y, U[0]);
			string fname = "p2_part3_" + to_string(t*del_t);
			PrintOut(fname + ".txt", U, n_x, n_y, del_x, del_y, 0);
			
		}
	}
	



	//clean up 
	for (size_t t = 0; t < 2; t++){
		for (size_t i = 0; i < n_x; i++){
			delete[] U[t][i];
		}
		delete[] U[t];
	}
	delete[]helper_container;
	delete[]b;

	
	
	
	return 0;


}
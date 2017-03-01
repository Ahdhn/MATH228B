//Solving problem 1 in homework three - MATH 228B - Winter 2017
//By Ahmed H. Mahmoud
//The code solves the problems, and produces raw data file for matlab to plot 

//#include "Plotter.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h> 
#include <algorithm>
#define _USE_MATH_DEFINES
#include<math.h>
#include <time.h>
using namespace std;



//Plotter _plot;//for plotting anywhere

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
inline void PrintOutTransAtPosition(double x,double y,double***temp,size_t n_x,size_t n_y,double del_x,double del_y,double del_t,size_t steady_time)
{
	//print out the transient solution (various times) for a certain position (x,y)
	size_t i,j,t;
	i=size_t(x/del_x);
	j=size_t(y/del_y);
	fstream file("Trans_Solution.txt",ios::out);
	file.precision(5);
	file<<"("<<x<<" "<<y<<")"<<endl;
	file<<"Time Temp"<<endl;

	for(t=0;t<steady_time;t++){
		file<<t*del_t<<" "<<temp[t][i][j]<<endl;		
	}
	file.close();

}

inline void ImposeBoundaryCondition(size_t n_x, double del_x, size_t n_y, double del_y, double***U, size_t time_size)
{
	//set the initial and boundary conditions 
	//at time>=0, U=0 @ 0 <= x <=X && y=0
	//at time>=0, U=0 @ 0 <  y < Y && x=0
	//at time>=0, U=0 @ 0 <= x <=X && y=Y
	//at time>=0, U=0 @ 0 <  y < Y && x=X
	
	//at time =0, T=exp(-10(x-0.3)^2 + (y-0.4)^2)) @ 0< x< X && 0<y<Y

	//1) initial conditions 	
	for (size_t i = 0; i<n_x; i++){
		for (size_t j = 0; j<n_y; j++){
			double x = double(i)*del_x;
			double y = double(j)*del_y;
			double xx = (x - 0.3)*(x - 0.3);
			double yy = (y - 0.4)*(y - 0.4);
			U[0][i][j] = exp(-10 * (xx + yy));
		}
	}

	//2) boundary conditions 
	//at time>=0, U=0 at 0 <= x <=X && y=0
	for (size_t t = 0; t < time_size; t++){
		for (size_t i = 0; i < n_x; i++){
			U[t][i][0] = 0;
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

inline void ADI(double***U, size_t n_x, size_t n_y, double del_x, double del_y, double del_t, size_t max_t, double*helper_container, double*b, double&my_time_mil_sec)
{
	//solving using ADI and reporting time in my_time_mil_sec
	//start with x-sweep then y-sweep
	//solving the tri-dia matrix using LU-decomposition method
	//http://www.math.buffalo.edu/~pitman/courses/mth437/na2/node2.html
	//at each sweep level, we only save two time steps
	//1) for x-sweep, we have n-time step at U[0] and save (n+0.5)-time step at U[1]
	//2) for y-sweep, we have (n+0.5)-time step at U[1] and save (n+1)-time step at U[0] 
	
	my_time_mil_sec = 0;

	if (n_x != n_y){
		std::cout << "Error (1) ADI(). n_x != n_y" << endl;				
		system("pause");
		exit(0);	
	}

	double r,alfa, beta, gamma; 
	
	r = (0.1*del_t) / (2.0*del_x*del_x);
	alfa = 1.0 + r;
	beta = 1.0 + 2.0*r;
	gamma = 1.0 - 2.0*r;
		

	clock_t t1, t2;
	t1 = clock();//get initial time 

	//fstream file("time_ittr.txt", ios::out);
	//file.precision(30);

	if (max_t == 641){
		PrintOut("0.txt", U, n_x, n_y, del_x, del_y, 0); 	
		//_plot.ColorIsoCountoring("0.ps", del_x, n_x, del_y, n_y, U[0]);
	}

	for (size_t t = 1; t < max_t; t++){ 
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

		double sum = 0;
		for (size_t i = 0; i < n_x; i++){
			for (size_t j = 0; j < n_y; j++){
				sum += U[0][i][j];
			}
		}

		if (max_t == 641 && (t % 50 == 0 )){ 
			string fname = to_string(t*del_t); 
			//_plot.ColorIsoCountoring(fname + ".ps", del_x, n_x, del_y, n_y, U[0]);			
			PrintOut(fname + ".txt", U, n_x, n_y, del_x, del_y, 0);			
		}
		//file << "time= " << t*del_t << " >>***<< sum_sol= " << sum << endl;
	}


	t2 = clock();//get final time 
	my_time_mil_sec += (double(t2) - double(t1)) / double(CLOCKS_PER_SEC * 1000.0);

	/*if (max_t == 641){
		PrintOut("1.txt", U, n_x, n_y, del_x, del_y, 0);
		_plot.ColorIsoCountoring("0.ps", del_x, n_x, del_y, n_y, U[0]);
	}*/	

}


int main (int argc, char *argv[])
{
	double***U;//solution at time step for all grid points 
	double X, Y, del_t, del_x, del_y, total_time;
	size_t n_x, n_y, max_t, n; 
	double*helper_container,*b;//helper array for tri-diag matrix solution 	

	//user input values	
	std::cout << "Solving problem 1 in homework three - MATH 228B - Winter 2017" << endl;

	total_time = 1.0;//total time 	
	X = 1.0;//length in x-direction
	Y = X;//length in y-direction
	
	del_x = 0.5;//initial grid spacing 
	del_t = 0.05; //initial time spacing 
	U = new double**[2];//we only need to store two time steps
	                    //either [n] and [n+(1/2)] or [n+(1/2)] and [n+1]

	//in this raw data file, we store the following 
	//del_x, del_y, del_t, time_mil_sec, total_sol
	fstream file("data.txt", ios::out);
	file.precision(30);
	file << "del_x  del_y del_t time_mil_sec total_sol" << endl;


	for (size_t ittr = 1; ittr < 10; ittr++){
		del_x /= 2.0;//each time halve the space step 
		del_t /= 2.0;//same thing for time step
		
		del_y = del_x;//same in x and y

		n_x = size_t(ceil(X / del_x) + 1); //gird num points 
		n_y = n_x;

		max_t = size_t(ceil(total_time / del_t) + 1);
		n = n_x*n_y; //total spacing points 
		

		//initlizing containers
		helper_container = new double[n - 1];
		b = new double[n - 1];

		for (size_t t = 0; t < 2; t++){
			U[t] = new double*[n_x];
			for (size_t i = 0; i < n_x; i++){
				U[t][i] = new double[n_y];
			}
		}

		ImposeBoundaryCondition(n_x, del_x, n_y, del_y, U, 2);

		double time_mil_sec;
		ADI(U, n_x, n_y, del_x, del_y, del_t, max_t, helper_container, b, time_mil_sec);

		//accumelate solution values for error calculation 
		double total_sol = 0;
		for (size_t i = 0; i < n_x; i++){
			for (size_t j = 0; j < n_y; j++){
				total_sol += U[0][i][j];
			}
		}
		
		//write to raw data file 
		file << del_x << "   " << del_y << "  " << del_t << "  " << time_mil_sec << "  " << total_sol << endl;
		
		//clean up for next iteration 
		for (size_t t = 0; t < 2; t++){
			for (size_t i = 0; i < n_x; i++){
				delete[] U[t][i];
			}
			delete[] U[t];
		}


		delete[]helper_container; 
		delete[]b;

	}
	
	
	return 0;


}
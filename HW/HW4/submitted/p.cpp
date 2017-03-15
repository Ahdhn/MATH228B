//Solving problem 1 in homework four - MATH 228B - Winter 2017
//By Ahmed H. Mahmoud
//Dependency: eigen library <http://eigen.tuxfamily.org/index.php?title=Main_Page>

//The code solves the problems, and produces raw data file for matlab to plot 
//The code carry on the refinement study using upwinding and Lax-Wendroff scheme with smooth and discontinious
//boundar condition over the 1d avection eqaution 

#include <iostream>
#include <fstream>
#include <string> 
#include <algorithm>
#define PI 3.14159265359
#include <Eigen\Dense>
using namespace Eigen;
size_t SCHEME_ID;
size_t IC_ID;

void PrintOut(std::string outputfilename, double*U, size_t n_x, double del_x)
{
	//print out raw data for matlab to draw 
	std::ofstream file(outputfilename, std::ios::out);
	file.precision(30);		
	file << n_x << std::endl;

	for (size_t i = 0; i < n_x; i++){ 			
		file << i*del_x << "	" << U[i] <<std::endl;
	}
	file.close();
	
}

double func(double x, bool smooth_condition)
{	
	
	if (smooth_condition){
		return cos(2.0*x*PI) + sin(6.0*x*PI) / 2.0;
	}
	else{
		if (abs(x - 0.5) < 0.25){
			return 1;
		}
		else{
			return 0;
		}		
	}
}
inline void ImposeBoundaryCondition(size_t n_x, double del_x, double*U, bool smooth_condition)
{
	//with smooth_condition, we use sin(x)
	//for discontinuous conditions, we use u = 1 @ (x-0.5)<0.25 
	//	                                       0 otherwise
	
	if (smooth_condition){
		for (size_t i = 0; i < n_x; i++){
			double x = double(i)*del_x;
			U[i] = func(x, smooth_condition);
		}
	}
	else{
		for (size_t i = 0; i < n_x; i++){
			double x = double(i)*del_x;
			U[i] = func(x, smooth_condition);
		
		}
	}
}

void Upwinding(double**U, double a, size_t n_x, double del_x, double del_t, size_t max_t)
{
	//Solve the 1d advection equation using the upwinding scheme with periodic boundary conditions 
	//a is the speed 
	double courant = (a*del_t) / del_x;	
	if (courant>1.0){
		std::cout << "Error (0) at Upwinding()::WARNING::Courant Number >1.0 and the scheme is unstable now!!!" << std::endl;
		system("pause");
	}
	double courant_min = 1 - courant;
	size_t n_plus = 1; 
	size_t n = 0;
	for (size_t t = 0; t < max_t; t++){
		for (size_t i = 0; i < n_x; i++){
			U[n_plus][i] = U[n][i] * courant_min;

			if (i == 0){ U[n_plus][i] += courant*U[n][n_x - 1]; }
			else{ U[n_plus][i] += courant*U[n][i- 1]; }
		}
		std::swap(n, n_plus);
	}
}
void LaxWendroff(double**U, double a, size_t n_x, double del_x, double del_t, size_t max_t)
{
	//Solve the 1d advection equation using the lax wendroff scheme with periodic boundary conditions 
	//a is the speed 
	double courant = a*del_t / del_x;
	
	if (courant>1.0){
		std::cout << "Error (0) at Upwinding()::WARNING::Courant Number >1.0 and the scheme is unstable now!!!" << std::endl;
		system("pause");
	}
	size_t n_plus = 1;
	size_t n = 0;

	double i_term = 1.0 - ((a*a*del_t*del_t) / (del_x*del_x)); 
	double i_mins_term = (( a*del_t) / (2.0*del_x)) + ((a*a*del_t*del_t) / (2.0*del_x*del_x)); 
	double i_plus_term = ((-a*del_t) / (2.0*del_x)) + ((a*a*del_t*del_t) / (2.0*del_x*del_x));

	for (size_t t = 0; t < max_t; t++){
		for (size_t i = 0; i < n_x; i++){
			U[n_plus][i] = U[n][i] * i_term;
			
			if (i == 0){ U[n_plus][i] += i_mins_term*U[n][n_x - 1]; }
			else{ U[n_plus][i] += i_mins_term*U[n][i - 1]; }

			if (i == n_x - 1){ U[n_plus][i] += i_plus_term*U[n][0]; }
			else{ U[n_plus][i] += i_plus_term*U[n][i + 1]; }

		}
		std::swap(n, n_plus);
	}
}

void SolveSys(int n, double A, double B, double C, double*b, double*x)
{
	//solving the square (n X n) tri-diagonla system with A as the main diagonal
	//C is the sub-diagona and B is the super diagonal
	//There is a C in the top right corner of the materix
	//and a B in the bottom left corner 
	//A B 0 0 0 C 
	//C A B 0 0 0 
	//0 C A B 0 0 
	//0 0 C A B 0 
	//0 0 0 C A B 
	//B 0 0 0 C A
	//b is the right hand side array
	//x is the unknown array
	
	//move data 
	MatrixXf A_mat(n, n);	
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			if (i == j){ 
				A_mat(i, j) = A; 
			}
			else if (i - j == -1){
				A_mat(i, j) = B; 
			}
			else if (j - i == -1){
				A_mat(i, j) = C; 
			}
			else{
				A_mat(i, j) = 0; 
			}
		}
	}
	A_mat(0, n - 1) = C;
	A_mat(n - 1, 0) = B; 

	VectorXf b_vec(n);
	VectorXf x_vec(n);
	for (int i = 0; i < n; i++){ 
		b_vec(i) = b[i]; 
		x_vec(i) = 0;
	}		

	x_vec = A_mat.colPivHouseholderQr().solve(b_vec);


	for (int i = 0; i < n; i++){		
		x[i] = x_vec(i);
	}	
}
inline void CrankNicolson(double**U, double a, size_t n_x, double del_x, double del_t, size_t max_t)
{
	//solving using Crank Nicolson's Scheme	
	//solving the tri-dia matix using the LU-decomposition method 
	//http://www.math.buffalo.edu/~pitman/courses/mth437/na2/node2.html
	//https://people.sc.fsu.edu/~jburkardt/classes/gateway_2014/lecture_week03.pdf

	//no need to declare the matix entities in vectors/arraies since they are all constants
	
	double courant = (a*del_t) / del_x;
	double A(1), B(courant*0.25), C(-courant*0.25);//tri-dia entities
	double* b = new double[n_x];
	size_t n_plus(1), n(0);	

	for (size_t t = 0; t < max_t; t++){ 
		//set up the RHS
		for (size_t i = 0; i < n_x; i++){ 			
			b[i] = U[n][i];
			if (i == 0){ b[i] += courant *0.25*U[n][n_x - 1]; }
			else{ b[i] += courant *0.25*U[n][i - 1]; }
			if (i == n_x - 1){ b[i] -= courant *0.25*U[n][0]; }
			else{ b[i] -= courant *0.25*U[n][i + 1]; }		
		}
		SolveSys(n_x, A, B, C, b, U[n_plus]);
		std::swap(n, n_plus);
	}
	
	free(b);
}

int main(int argc, char *argv[])
{
	
	std::cout << "Solving problem 1 in homework four - MATH 228B - Winter 2017" << std::endl;
	IC_ID = 5;
	while (IC_ID != 0 && IC_ID != 1){
		std::cout << "\nEnter 1 for smooth initial conditions or 0 for discontinuous initial conditions" << std::endl;
		std::cin >> IC_ID;
		if (IC_ID != 0 && IC_ID != 1){
			std::cout << "Invalid Entry!!!" << std::endl;
		}
	}
	
	SCHEME_ID = 5;
	while (SCHEME_ID != 1 && SCHEME_ID != 2 && SCHEME_ID != 3){
		std::cout << "\nEnter 1 for Upwinding scheme, 2 for Lax-Wendroff scheme, or 3 for Crank-Nicolson-like scheme" << std::endl;
		std::cin >> SCHEME_ID;
		if (SCHEME_ID != 1 && SCHEME_ID != 2 && SCHEME_ID != 3){
			std::cout << "Invalid Entry!!!" << std::endl;
		}
	}

	double**U;//solution at time step for all grid points 
	double a = 1.0;//speed 

	double total_time = 1.0;//total time 	
	double X = 1.0;//length in x-direction
	
	
	double del_x = 0.5;//initial grid spacing 	
	U = new double*[2];//we only need to store two time steps
	                    //either [n] and [n+(1/2)] or [n+(1/2)] and [n+1]


	//in this raw data file, we store the following 
	//del_x, del_t, time_mil_sec, total_sol
	std::string fname;
	if (SCHEME_ID == 1){ fname = "UP/Upwinding_data.txt"; }
	else if (SCHEME_ID == 2){ fname = "LW/LW_data.txt"; }
	else if (SCHEME_ID == 3){ fname = "CN/CN_data.txt"; }
	std::fstream file(fname, std::ios::out);
	file.precision(30);
	file << "del_t		n_t		del_x		n_x		max_err		1_norm		2_norm" << std::endl;


	for (size_t ittr = 1; ittr < 15; ittr++){
		del_x /= 2.0;//each time halve the space step 

		if (SCHEME_ID == 3){ del_x = 1.0 / 1000.0; }//jump straight ahead to fine mesh for CN because it is done once 

		double del_t = 0.9*a*del_x;

		size_t n_x = size_t(ceil(X / del_x) + 1); //gird num points 
		size_t max_t = size_t(ceil(total_time / del_t) + 1);//number of time steps 
		
		//initlizing containers
		U[0] = new double[n_x];
		U[1] = new double[n_x];
		
		ImposeBoundaryCondition(n_x, del_x, U[0], IC_ID==1);

		//print out the ground truth 
		/*std::string  fname_ground = std::to_string(ittr) + "_GroundTruth_X" + std::to_string(del_x) + "_T" + std::to_string(del_t) + ".txt";
		if (SCHEME_ID == 1){ fname_ground = "UP/" + std::to_string(ittr) + "_GroundTruth_X" + std::to_string(del_x) + "_T" + std::to_string(del_t) + ".txt"; }
		else if (SCHEME_ID == 2){ fname_ground = "LW/" + std::to_string(ittr) + "_GroundTruth_X" + std::to_string(del_x) + "_T" + std::to_string(del_t) + ".txt"; }
		else if (SCHEME_ID == 3){ fname_ground = "CN/" + std::to_string(ittr) + "_GroundTruth_X" + std::to_string(del_x) + "_T" + std::to_string(del_t) + ".txt"; }
		PrintOut(fname_ground, U[0], n_x, del_x);*/


		if (SCHEME_ID == 1){
			Upwinding(U, a, n_x, del_x, del_t, max_t);
		}
		else if (SCHEME_ID == 2){
			LaxWendroff(U, a, n_x, del_x, del_t, max_t);
		}
		else if (SCHEME_ID == 3){
			CrankNicolson(U, a, n_x, del_x, del_t, max_t); 
		}
		else{
			std::cout << "Error:: Invalid SHCEME_ID " << std::endl;
			system("pause");
		}


		size_t n_plus = (max_t % 2 == 0) ? 0 : 1; //where the last update has been placed 
		//calc error (max, 1-norm, 2-norm)
		double max_err(-1.0),sum_err(0), sum_err_sq(0);

		for (size_t i = 0; i < n_x; i++){	
			double x = double(i)*del_x;
			double exact = func(x, IC_ID == 1);
			double myErr = abs((U[n_plus][i] - exact)); 
			

			max_err = std::max(max_err, myErr);
			sum_err += myErr;
			sum_err_sq += myErr*myErr;
		}

		sum_err *= del_x;//1-norm
		sum_err_sq = sqrt(del_x*sum_err_sq); 
		//write to raw data file 
		file << del_t << "		" << max_t << "		" << del_x << "		" << n_x << "		" << max_err << "		" << sum_err << "		" << sum_err_sq << std::endl; 
		
		/*std::string fname_matlab;
		if (SCHEME_ID == 1){ fname_matlab = "UP/" + std::to_string(ittr) + "_Upwinding_X" + std::to_string(del_x) + "_T" + std::to_string(del_t) + ".txt"; }
		else if (SCHEME_ID == 2){ fname_matlab = "LW/" + std::to_string(ittr) + "_LW_X" + std::to_string(del_x) + "_T" + std::to_string(del_t) + ".txt"; }
		else if (SCHEME_ID == 3){ fname_matlab ="CN/"+ std::to_string(ittr) + "_CN_X" + std::to_string(del_x) + "_T" + std::to_string(del_t) + ".txt"; }
		PrintOut(fname_matlab, U[n_plus], n_x, del_x);*/


		//clean up for next iteration 
		delete[] U[0];
		delete[] U[1];

		if (SCHEME_ID == 3){ break; }//only once for CN
	}
	
	delete[] U;
	return 0;


}
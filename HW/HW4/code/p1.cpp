//Solving problem 1 in homework four - MATH 228B - Winter 2017
//By Ahmed H. Mahmoud
//The code solves the problems, and produces raw data file for matlab to plot 
//The code carry on the refinement study using upwinding and Lax-Wendroff scheme with smooth and discontinious
//boundar condition over the 1d avection eqaution 

#include <iostream>
#include <fstream>
#include <string> 
#include <algorithm>
#define PI 3.14159265359


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

double func(double x)
{
	
	
	//return sin(x*(2.0*PI));	
	return cos(2.0*x*PI) + sin(6.0*x*PI) / 2.0; 
}
inline void ImposeBoundaryCondition(size_t n_x, double del_x, double*U, bool smooth_condition)
{
	//with smooth_condition, we use sin(x)
	//for discontinuous conditions, we use u = 1 @ (x-0.5)<0.25 
	//	                                       0 otherwise
	
	if (smooth_condition){
		for (size_t i = 0; i < n_x; i++){
			double x = double(i)*del_x;
			U[i] = func(x);
		}
	}
	else{
		for (size_t i = 0; i < n_x; i++){
			double x = double(i)*del_x;
			if (x - 0.5 < 0.25){
				U[i] = 1;
			}
			else{
				U[i] = 0;
			}
		
		}
	}
}

void Upwinding(double**U, double a, size_t n_x, double del_x, double del_t, size_t max_t)
{
	//Solve the 1d advection equation using the upwinding scheme with periodic boundary conditions 
	//a is the speed 
	double courant = a*del_t / del_x;
	double courant_min = 1 - courant ;
	if (courant>1.0){
		std::cout << "Error (0) at Upwinding()::WARNING::Courant Number >1.0 and the scheme is unstable now!!!" << std::endl;
		system("pause");
	}
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
	double i_mins_term = ((a*del_t) / (2.0*del_x)) + ((a*a*del_t*del_t) / (2.0*del_x*del_x)); 
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
void CN_like(double**U, double a, size_t n_x, double del_x, double del_t, size_t max_t)
{
	//Solve the 1d advection equation using the the CN-like scheme with periodic boundary conditions 
	//a is the speed 
		
	size_t n_plus = 1;
	size_t n = 0;

	double i_term = 1.0 - ((a*a*del_t*del_t) / (del_x*del_x));
	double i_mins_term = ((a*del_t) / (2.0*del_x)) + ((a*a*del_t*del_t) / (2.0*del_x*del_x));
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
int main (int argc, char *argv[])
{
	std::cout << "Solving problem 1 in homework four - MATH 228B - Winter 2017" << std::endl;
	IC_ID = 5;
	while (IC_ID != 0 && IC_ID != 1){
		std::cout << "\nEnter 1 for smooth initial conditions (sin(x)) or 0 for discontinuous initial conditions" << std::endl;
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
	
	
	double del_x = 1.0;//initial grid spacing 	
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
		double del_t = 0.9*a*del_x;

		size_t n_x = size_t(ceil(X / del_x) + 1); //gird num points 
		size_t max_t = size_t(ceil(total_time / del_t) + 1);//number of time steps 
		
		//initlizing containers
		U[0] = new double[n_x];
		U[1] = new double[n_x];
		
		ImposeBoundaryCondition(n_x, del_x, U[0], IC_ID==1);

		//print out the ground truth 
		std::string  fname_ground = std::to_string(ittr) + "_GroundTruth_X" + std::to_string(del_x) + "_T" + std::to_string(del_t) + ".txt";
		if (SCHEME_ID == 1){ fname_ground = "UP/" + std::to_string(ittr) + "_GroundTruth_X" + std::to_string(del_x) + "_T" + std::to_string(del_t) + ".txt"; }
		else if (SCHEME_ID == 2){ fname_ground = "LW/" + std::to_string(ittr) + "_GroundTruth_X" + std::to_string(del_x) + "_T" + std::to_string(del_t) + ".txt"; }
		else if (SCHEME_ID == 3){ fname_ground = "CN/" + std::to_string(ittr) + "_GroundTruth_X" + std::to_string(del_x) + "_T" + std::to_string(del_t) + ".txt"; }
		PrintOut(fname_ground, U[0], n_x, del_x);


		if (SCHEME_ID == 1){
			Upwinding(U, 1.0, n_x, del_x, del_t, max_t);
		}
		else if (SCHEME_ID == 2){
			LaxWendroff(U, 1.0, n_x, del_x, del_t, max_t);
		}
		else if (SCHEME_ID == 3){
			CN_like(U, 1.0, n_x, del_x, del_t, max_t);
		}
		else{
			std::cout << "Error:: Invalid SHCEME_ID " << std::endl;
			system("pause");
		}


		size_t n_plus = (max_t % 2 == 0) ? 0 : 1; //where the last update has been placed 
		//calc error (max, 1-norm, 2-norm)
		double max_err(0),sum_err(0), sum_err_sq(0);

		for (size_t i = 0; i < n_x; i++){	
			double x = double(i)*del_x;
			double myErr = abs(U[n_plus][i] - func(x));

			max_err = std::max(max_err, myErr);
			sum_err += myErr;
			sum_err_sq += myErr*myErr;
		}

		
		//write to raw data file 
		file << del_t << "		" << max_t << "		" << del_x << "		" << n_x << "		" << max_err << "		" << sum_err / n_x << "		" << sqrt(sum_err_sq) / n_x << std::endl; 
		
		std::string fname_matlab;
		if (SCHEME_ID == 1){ fname_matlab = "UP/" + std::to_string(ittr) + "_Upwinding_X" + std::to_string(del_x) + "_T" + std::to_string(del_t) + ".txt"; }
		else if (SCHEME_ID == 2){ fname_matlab = "LW/" + std::to_string(ittr) + "_LW_X" + std::to_string(del_x) + "_T" + std::to_string(del_t) + ".txt"; }
		else if (SCHEME_ID == 3){ fname_matlab ="CN/"+ std::to_string(ittr) + "_CN_X" + std::to_string(del_x) + "_T" + std::to_string(del_t) + ".txt"; }
		PrintOut(fname_matlab, U[n_plus], n_x, del_x);


		//clean up for next iteration 
		delete[] U[0];
		delete[] U[1];
	}
	
	delete[] U;
	return 0;


}
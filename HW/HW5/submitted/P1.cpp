//Solving problem 1 in homework five - MATH 228B - Winter 2017
//By Ahmed H. Mahmoud

//The code solves the one spatial dimesion of the linearized equations of acoustics 
//(sound waves) using Lax-Wendroff on [0,1] using cell-centered grid 

#include <iostream>
#include <fstream>
#include <string> 
#include <algorithm>
#define PI 3.14159265359

int IC_ID;

void PrintOut(std::string outputfilename, double*U, size_t n_x, double del_x)
{
	//print out raw data for matlab to draw 
	std::ofstream file(outputfilename, std::ios::out);
	file.precision(30);
	//file << n_x << std::endl;

	for (size_t i = 0; i < n_x; i++){
		file << i*del_x << "	" << U[i] << std::endl;
	}
	file.close();

}
void WriteData(double*P, double*U,size_t n_x,double del_x, double time)
{
	//Write data out 
	std::string fname_matlab_P, fname_matlab_U;
	std::string ic_str;
	if (IC_ID == 1){
		ic_str = "_sin_";
	}
	else if (IC_ID == 2){
		ic_str = "_tri_";
	}
	else if (IC_ID == 3){
		ic_str = "_sin_tri_";
	}
	else if (IC_ID == 4){
		ic_str = "_sin_expo_";
	}
	fname_matlab_P = "P_X" + ic_str + std::to_string(time) + ".txt";
	fname_matlab_U = "U_X" + ic_str + std::to_string(time) + ".txt";
	PrintOut(fname_matlab_P, P, n_x, del_x);
	PrintOut(fname_matlab_U, U, n_x, del_x);
}
void ImposeBoundaryCondition(double*U, double*P, double K, double Rho, double Courant, size_t n_x, double del_x)
{
	//Initial conditions are 
	//U is always zero
	//IC_ID ==1 --> sin(10 PI x) @ 0.4 <= x < 0.6
	//IC_ID ==2 --> sin(40 PI x) @ 0.4 <= x < 0.6
	//IC_ID ==3 --> 1 - x        @ x < 1
	//              0            otherwise 
	//IC_ID ==4 --> expo

	for (size_t i = 0; i < n_x; i++){
		U[i] = 0;
		double x = double(i)*del_x;
		if (IC_ID == 1){
			if (x >= 0.4 && x <= 0.6){				
				P[i] = sin(10 * PI*x); 
			}
			else{
				P[i] = 0;
			}
		}
		else if (IC_ID == 3){
			if (x >= 0.2 && x <= 0.4){
				P[i] = sin(10 * PI*x);
			}
			else if (x >= 0.6 && x <= 0.8){
				P[i] = 10.0*(0.1 - abs(x - 0.7));
			}
			else{
				P[i] = 0;
			}
		}
		else if (IC_ID == 2){
			if (x >= 0.2 && x <= 0.8){
				P[i] = (10.0 / 3.0)*(0.3 - abs(x - 0.5)); 
			}
			else{
				P[i] = 0;
			}
		}
		else if (IC_ID == 4){
			P[i] = exp(-5*x)*cos(10 * PI*x);
			/*if (x >= 0.4 && x <= 0.6){
				P[i] = 0;
			}
			else{
				P[i] = 0;
			}*/
		}
		else{
			std::cout << "Error::ImposeBoundaryCondition" << std::endl;
			system("pause");
		}
	}

	//ghost points 
	U[0] = -U[1]; 
	P[0] = P[1];

	U[n_x - 1] = 0.5*(U[n_x - 2] + (P[n_x - 2] / sqrt(K*Rho))); 
	P[n_x - 1] = 0.5*(P[n_x - 2] + (U[n_x - 2] * sqrt(K*Rho))); 
	
}

void LaxWendroff(double ** U, double**P, double K, double Rho, double Courant, size_t n_x, double del_x, double del_t, size_t max_t)
{
	//solve the 1d linearized equation of sound waves using lax wendroff method 
	if (Courant>1.0){
		std::cout << "Error (0) at LaxWendroff()::WARNING::Courant Number >1.0 and the scheme is unstable now!!!" << std::endl;
		system("pause");
	}

	size_t n_plus = 1;
	size_t n = 0;
	double inv_rho = 1 / Rho;
	double i_term_p = (-K*K *del_t*del_t) / (del_x*del_x);
	double i_term_u = (-inv_rho*inv_rho*del_t*del_t) / (del_x*del_x); 
	double i_mins_p = ((K*K*del_t*del_t) / (2.0*del_x*del_x)) + ((K*del_t) / (2.0*del_x));
	double i_mins_u = ((inv_rho*inv_rho*del_t*del_t) / (2.0*del_x*del_x)) + ( (inv_rho*del_t) / (2.0*del_x));
	double i_plus_p = ((K*K*del_t*del_t) / (2.0*del_x*del_x)) - ((K*del_t) / (2.0*del_x)); 
	double i_plus_u = ((inv_rho*inv_rho*del_t*del_t) / (2.0*del_x*del_x)) - ((inv_rho*del_t) / (2.0*del_x));
	
	double term1_p = (-del_t*K) / (2.0*del_x);
	double term1_u = (-del_t*inv_rho) / (2.0*del_x);
	double term2 = (del_t*del_t*K) / (2.0*del_x*del_x*Rho);
	

	for (size_t t = 0; t < max_t; t++){
		if (t % 50 == 0){
			WriteData(P[n], U[n], n_x, del_x, t*del_t);
		}

		for (size_t i = 1; i < n_x - 1; i++){//the first (i=0) and last (i=n_x-1) are ghost points
			//P[n_plus][i] = P[n][i] + i_term_p*U[n][i] + i_mins_p*U[n][i - 1] + i_plus_p*U[n][i + 1]; 
			//U[n_plus][i] = U[n][i] + i_term_u*P[n][i] + i_mins_u*P[n][i - 1] + i_plus_u*P[n][i + 1];

			P[n_plus][i] = P[n][i] + term1_p*(U[n][i + 1] - U[n][i - 1]) + term2*(P[n][i + 1] - 2.0*P[n][i] + P[n][i - 1]);
			U[n_plus][i] = U[n][i] + term1_u*(P[n][i + 1] - P[n][i - 1]) + term2*(U[n][i + 1] - 2.0*U[n][i] + U[n][i - 1]);
			
		}
		U[n_plus][0] = -U[n_plus][1];
		P[n_plus][0] = P[n_plus][1];

		U[n_plus][n_x - 1] = 0.5*(U[n_plus][n_x - 2] + (P[n_plus][n_x - 2] / sqrt(K*Rho)));
		P[n_plus][n_x - 1] = 0.5*(P[n_plus][n_x - 2] + (U[n_plus][n_x - 2] * sqrt(K*Rho)));

		

		std::swap(n, n_plus);
	}

}

int main(int argc, char *argv[])
{

	std::cout << "Solving problem 1 in homework five - MATH 228B - Winter 2017" << std::endl;
	IC_ID = 5;
	while (IC_ID != 1 && IC_ID != 2 && IC_ID != 3 && IC_ID != 4){ 
		std::cout << "\nPlease enter: \n1 for sin wave \n2 for triangular function \n3 for sin and triangular function \n4 for decaying sinusoidal exponential function" << std::endl;
		std::cin >> IC_ID;
		if (IC_ID != 1 && IC_ID != 2 && IC_ID != 3 && IC_ID != 4){
			std::cout << "Invalid Entry!!!" << std::endl;
		}
	}


	double**U,**P;//solution at time step for all grid points 
	double K(0.3), Rho(1.0), Courant(0.8);

	double total_time = 5.0;//total time 	
	double X = 1.0;//length in x-direction


	double del_x = 0.005;
	double del_t = (Courant*del_x)/* / (sqrt(K / Rho))*/;

	//double del_t = 1.0 / 320.0;
	//double del_x = Courant *del_t;

	U = new double*[2];//we only need to store two time steps
	P = new double*[2];//either [n] and [n+(1/2)] or [n+(1/2)] and [n+1]

	size_t n_x = size_t(ceil(X / del_x) + 3); //gird num points + ghost points 
	size_t max_t = size_t(ceil(total_time / del_t) + 1);//number of time steps 

	//initlizing containers
	U[0] = new double[n_x];
	U[1] = new double[n_x];
	P[0] = new double[n_x];
	P[1] = new double[n_x];

	ImposeBoundaryCondition(U[0], P[0], K, Rho, Courant, n_x, del_x); 
	//std::string fname_matlab;
	//fname_matlab = "init.txt";
	//PrintOut(fname_matlab, P[0], n_x, del_x);

	LaxWendroff(U,P,K,Rho,Courant, n_x, del_x, del_t, max_t);
	
	size_t n_plus = (max_t % 2 == 0) ? 0 : 1;
	WriteData(P[n_plus], U[n_plus], n_x, del_x, total_time);




	return 0;
}
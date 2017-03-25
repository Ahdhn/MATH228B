//Solving problem 2 in homework five - MATH 228B - Winter 2017
//By Ahmed H. Mahmoud
 

#include <iostream>
#include <fstream>
#include <string> 
#include <algorithm>
#define PI 3.14159265359

int IC_ID;
int LIMITER_ID;

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

void ImposeBoundaryCondition(double*U, size_t n_x, double del_x)
{
	//Initial conditions are 
	//U is always zero
	//IC_ID ==1 --> wave packet u(x,0) = cos(16PI*x) exp(-50*(x-0.5)^2)
	//IC_ID ==2 --> smooth, low frequency u(x,0) = sin(2PI*x)sin(4*PI*x)
	//IC_ID ==3 --> step function  1    if |x-0.5|<0.25
	//                             0    Otherwise 
	

	for (size_t i = 0; i < n_x; i++){		
		double x = double(i)*del_x;
		if (IC_ID == 1){					
			U[i] = cos(16 * PI*x)*exp(-50 * (x - 0.5)*(x - 0.5));			
		}
		else if (IC_ID == 2){
			U[i] = sin(2 * PI*x)*sin(4 * PI*x);
		}
		else if (IC_ID == 3){
			if (abs(x - 0.5) < 0.25){
				U[i] = 1.0;
			}
			else{
				U[i] = 0.0;
			}		
		}		
		else{
			std::cout << "Error::ImposeBoundaryCondition" << std::endl;
			system("pause");
		}
	}


}
inline double Limiter(int id, double theta)
{
	//return the flux limiter phi(theta)
	if (id == 1){
		//Upwinding
		return 0;
	}
	else if (id == 2){ 
		//LW
		return 1;
	}
	else if (id == 3){
		//BW
		return theta;
	}
	else if (id == 4){
		//minmod
		return std::max(0.0, std::min(1.0, theta));
	}
	else if (id == 5){
		//superbee
		double inter = std::max(0.0, std::min(1.0, 2.0*theta));
		return  std::max(inter, std::min(2.0, theta));
	}
	else if (id == 6){
		//MC
		double inter = std::min((1.0 + theta) / 2.0, 2.0);
		return std::max(0.0, std::min(inter, 2.0*theta));
	}
	else if (id == 7){
		//van Leer
		return (theta + abs(theta)) / (1.0 + abs(theta));
	}
	else {
		std::cout << "Error::Limiter()" << std::endl;
		system("pause");
	}
}
double Flux(double*U, size_t j, double a, double del_t, double del_x, size_t n_x)
{
	//get the flux 
	size_t j_mins = (j == 0) ? n_x - 1 : j - 1;
	size_t j_plus = (j == n_x-1) ? 0 : j + 1;

	double diff = U[j] - U[j_mins]; 

	double f_up = (a >= 0) ? a*U[j_mins] : a*U[j];
	size_t J_up = (a >= 0) ? j_mins : j_plus;
	
	size_t J_up_minus = (J_up == 0) ? n_x - 1 : J_up - 1;
	double theta = (U[J_up] - U[J_up_minus]) / (U[j] - U[j_mins]);
	
	double delta_j_minus = diff*Limiter(LIMITER_ID, theta);

	return f_up + 0.5*abs(a)*delta_j_minus*(1 - abs(a*del_t / del_x));
}
void Solver(double ** U, double a, double Courant, double limiter_id, size_t n_x, double del_x, double del_t, size_t max_t)
{
	size_t n_plus = 1;
	size_t n = 0;
	double ratio = del_t / del_x;

	for (size_t t = 0; t < max_t; t++){
		for (size_t i = 0; i < n_x; i++){//the first (i=0) and last (i=n_x-1) are ghost points
			double F_j_mius = Flux(U[n], i, a, del_t, del_x,n_x); 
			double F_j_plus = Flux(U[n], (i==n_x-1) ? 0 :  i + 1, a, del_t, del_x, n_x); //sending j+1 because this function return j-0.5, thus j+1-0.5 = j+0.5
			U[n_plus][i] = U[n][i] - ratio*(F_j_plus - F_j_mius);
			if (U[n_plus][i]<-1.0){
				size_t dfgdf = 54545;
			}
		}		
		std::swap(n, n_plus);
	}

}
int main(int argc, char *argv[])
{

	std::cout << "Solving problem 2 in homework five - MATH 228B - Winter 2017" << std::endl;
	IC_ID = 5;
	while (IC_ID != 1 && IC_ID != 2 && IC_ID != 3){
		std::cout << "\nFor initial conditions, please enter: \n1 for Wave packet \n2 for smooth, low frequency \n3 for step function function" << std::endl;
		std::cin >> IC_ID;
		if (IC_ID != 1 && IC_ID != 2 && IC_ID != 3){
			std::cout << "Invalid Entry!!!" << std::endl;
		}
	}

	LIMITER_ID == 100;
	while (LIMITER_ID != 1 && LIMITER_ID != 2 && LIMITER_ID != 3 && LIMITER_ID != 4 && LIMITER_ID != 5 && LIMITER_ID != 6 && LIMITER_ID != 7){
		std::cout << "\nFor limiter choice, please enter: \n1 for Upwinding \n2 for Lax-Wendroff \n3 for Beam-Warming \n4 for minmod, \n5 for Superbee, \n6 for MC, \n7 for van Leer" << std::endl;
		std::cin >> LIMITER_ID;
		if (LIMITER_ID != 1 && LIMITER_ID != 2 && LIMITER_ID != 3 && LIMITER_ID != 4 && LIMITER_ID != 5 && LIMITER_ID != 6 && LIMITER_ID != 7){
			std::cout << "Invalid Entry!!!" << std::endl;
		}
	}

	double**U;//solution at time step for all grid points 
	double a(1.0), Courant(0.9);

	double total_time = 5.0;//total time 	
	double X = 1.0;//length in x-direction


	double del_x = 0.005;
	double del_t = (Courant*del_x) / a;


	U = new double*[2];//we only need to store two time steps

	size_t n_x = size_t(ceil(X / del_x)); //gird num points + ghost points 
	size_t max_t = size_t(ceil(total_time / del_t) );//number of time steps 

	//initlizing containers
	U[0] = new double[n_x];
	U[1] = new double[n_x];

	ImposeBoundaryCondition(U[0], n_x, del_x);
	std::string fname_matlab;
	fname_matlab = "init.txt";
	PrintOut(fname_matlab, U[0], n_x, del_x);

	
	Solver(U, a, Courant, LIMITER_ID, n_x, del_x, del_t, max_t);


	//Write data out 
	std::string fname_matlab_U;
	std::string ic_str;
	if (IC_ID == 1){
		ic_str = "_Wave_";
	}
	else if (IC_ID == 2){
		ic_str = "_Low_freq_";
	}
	else if (IC_ID == 3){
		ic_str = "_Step_";
	}
	
	std::string lim_str;
	if (LIMITER_ID == 1){
		lim_str = "Upwinding_";
	}
	else if (LIMITER_ID == 2){
		lim_str = "LW_";
	}
	else if (LIMITER_ID == 3){
		lim_str = "BW_";
	}
	else if (LIMITER_ID == 4){
		lim_str = "Minmod_";
	}
	else if (LIMITER_ID == 5){
		lim_str = "Superbee_";
	}
	else if (LIMITER_ID == 6){
		lim_str = "MC_";
	}
	else if (LIMITER_ID == 7){
		lim_str = "vanLeer_";
	}


	fname_matlab_U = "U" + ic_str + lim_str + ".txt";
	size_t n_plus = (max_t % 2 == 0) ? 0 : 1;	
	PrintOut(fname_matlab_U, U[n_plus], n_x, del_x);




	return 0;
}
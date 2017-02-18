//Solving problem 2 and 3 in homework one - MATH 228B - Winter 2017
//By Ahmed H. Mahmoud
//The code solves the problems, and produces raw data file for matlab to plot 

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h> 
#define _USE_MATH_DEFINES
#include<math.h>

//#include "Plotter.h"

using namespace std;

size_t _PROBLEM;//probelm flag

inline void PrintOut(string outputfilename,double**U,size_t n_t,size_t n_x,double del_t,double del_x)
{
	//print out the U array (U at all time in different sections) in a txt file
	//to be used for plotting

	ofstream file(outputfilename, ios::out);
	file.precision(10);
	size_t t, i;
	double x, time;
	file << n_x << " " << n_t << endl;
	for (t = 0; t<n_t; t++){
		time = t*del_t;
		for (i = 0; i<n_x; i++){
			x = i*del_x;
			file << U[t][i] << " ";
		}
		file << endl;
	}
}

inline void Decomp (size_t n_x, double alfa, double gamma, double beta, double*l, double*Mu)
{
	//decompose the tri-dia matrix into two matrices
	//http://www.math.buffalo.edu/~pitman/courses/mth437/na2/node3.html
	//http://www.math.buffalo.edu/~pitman/courses/mth437/na2/node4.html
	
	l[1]=alfa;
	for (size_t i = 2; i<n_x - 1; i++){
		Mu[i - 1] = gamma / l[i - 1];		
		l[i] = alfa - beta*Mu[i - 1];		
	}
}
inline void Solv(size_t n_x,double*b, double beta, double*l, double*Mu, double*z, double*&X)
{
	//solve the decomposed matirces 
	
	z[1] = b[1] / l[1];
	for (size_t i = 2; i<n_x - 1; i++){
		z[i] = (b[i] - beta*z[i - 1]) / l[i];
	}

	X[n_x - 2] = z[n_x - 2];
	for (size_t i = n_x - 3; i >= 1; i--){
		X[i] = z[i] - Mu[i] * X[i + 1];
	}
}

inline void ImposeBoundaryCondition(size_t n_t, size_t n_x, double**U, double del_t, double del_x)
{
	//set the initial and boundary conditions according to the problem 
	//for _PROBLEM == 2
	//at time = 0, U=0 for 0<x<X
	//at time > 0, U=0 at x=0 and x=X
	//for _PROBLEM == 3 || _PROBLEM == 4
	//at time = 0, U=1 for 0<x<0.5, U=0 for 0.5<x<X
	//at time > 0, U=1 at x=0 and U=0 at x=X

	size_t i;
	if (_PROBLEM == 3 || _PROBLEM == 4){
		for (i = 0; i < n_t; i++){ 
			U[i][0] = 1;
			U[i][n_x - 1] = 0;
		}
		for (i = 1; i < n_x; i++){ 
			if (i*del_x < 0.5){
				U[0][i] = 1;
			}
			else{
				U[0][i] = 0;
			}
		}
	}
	else if (_PROBLEM == 2){
		for (i = 0; i<n_t; i++){
			U[i][0] = 0;			
			U[i][n_x-1] = 0;
		}
		for (i = 1; i<n_x; i++){
			U[0][i] = 0;
		}
	}
}


inline void CrankNicolson(double**U,size_t n_t,size_t n_x,double alpha,double del_t,double del_x,
	                      double*l, double*Mu, double*z, double*b)
{
	//solving using Crank Nicolson's Scheme	
	//solving the tri-dia matix using the LU-decomposition method 
	//http://www.math.buffalo.edu/~pitman/courses/mth437/na2/node2.html
	//https://people.sc.fsu.edu/~jburkardt/classes/gateway_2014/lecture_week03.pdf

	//no need to declare the matix entities in vectors/arraies since they are all constants

	double cof = 0.5*(alpha*del_t) / (del_x*del_x); //coefficient 
	double alfa(2 * cof + 1), beta(-cof), gamma(-cof);//tri-dia entities
	
		
	for (size_t t = 1; t < n_t; t++){ 
		//matrix setup for this time step
		for (size_t V = 1; V<n_x - 1; V++){
			b[V] = (1 - 2 * cof)*U[t - 1][V] + cof*(U[t - 1][V + 1] + U[t - 1][V - 1]);
			if (_PROBLEM == 2){
				b[V] -= 0.5*del_t*(2.0 - exp(-(double(t) - 1)*del_t) - exp(-double(t)*del_t));
			}
		}
		b[1] -= gamma*U[t][0];
		b[n_x - 2] -= beta*U[t][n_x - 1];

		//solve
		Decomp(n_x, alfa, gamma, beta, l, Mu);
		Solv(n_x, b, beta, l, Mu, z, U[t]);
	}

}

inline void BDF2(double**U, size_t n_t, size_t n_x, double alpha, double del_t, double del_x,
	             double*l, double*Mu, double*z, double*b)
{
	//solving using BDF2's Scheme	
	//using BDF1 to slove the first time step, then carry on the solution on BDF2
	//solving the tri-dia matix using the LU-decomposition method 
	
	//no need to declare the matix entities in vectors/arraies since they are all constants

	double cof = (alpha*del_t) / (del_x*del_x); //coefficient 
	double alfa2(4 * cof + 3), beta2(-2.0*cof), gamma2(-2.0*cof);//tri-dia entities for BFD2
	double alfa1(2 * cof + 1), beta1(-cof), gamma1(-cof);//tri-dia entities for BFD1
	
	//BDF1 for the first step in time 
	for (size_t t = 1; t < 2; t++){//just one step
		for (size_t V = 1; V < n_x - 1; V++){
			b[V] = U[t - 1][V];
		}
		b[1] -= gamma1*U[t][0];
		b[n_x - 2] -= beta1*U[t][n_x - 1];

		//solve
		Decomp(n_x, alfa1, gamma1, beta1, l, Mu);
		Solv(n_x, b, beta1, l, Mu, z, U[t]);
	}


	//BDF2 for the rest time steps 
	for (size_t t = 2; t < n_t; t++){ 
		//matrix setup for this time step
		for (size_t V = 1; V<n_x - 1; V++){
			b[V] = 4*U[t - 1][V] - U[t - 2][V];			
		}
		b[1] -= gamma2*U[t][0];
		b[n_x - 2] -= beta2*U[t][n_x - 1];

		//solve
		Decomp(n_x, alfa2, gamma2, beta2, l, Mu);
		Solv(n_x, b, beta2, l, Mu, z, U[t]);
	}

}

int main (int argc, char *argv[])
{
	double X;//total spcae
	double del_t,del_x,time;//time step and speca step, total time
	size_t n_x,n_t; //number of elements in x and t direction
	double**U;//solution for n times to all grid points
	double*l, *Mu, *z, *b;//helper array for tri-diag matrix solution 	

	//user input values
	cout<<"Solving problem 1 and 2 in homework one - MATH 228B - Winter 2017"<<endl;	

	_PROBLEM = 0;

	cout<<"\n******** PROBLEM ********"<<endl;	
	while(_PROBLEM!=2 && _PROBLEM!=3){		
		cout<<"\nPlease enter 2 for problem no.2 or 3 for problem no.3: "<<endl;		
		cin>>_PROBLEM;
		if(_PROBLEM!=2 && _PROBLEM!=3){
			cout<<"Invalid Entry!!!"<<endl;
		}
		else if (_PROBLEM == 3){
			cout << "\nPlease enter 1 for Cranck-Nicolson or 2 for BDF2" << endl;
			cin >> _PROBLEM;
			if (_PROBLEM != 1 && _PROBLEM != 2){
				cout << "Invalid Entry!!!" << endl;
			}
			else{
				if (_PROBLEM == 1){ _PROBLEM = 3; }
				else{ _PROBLEM = 4; }
				break;
			}
		}
	}

	cout<<"\n******** GEOMETRY & BOUNDARY CONDITIONS ********"<<endl;
	time = 1;	//total time 	
	X=1; 	//total space		
	del_t = 1.0 / 10.0;//time step 	
	del_x = 1.0 / 50.0; //space step 
	
	n_x = 1 + size_t(X / del_x);
	n_t = 1 + size_t(time / del_t);

	//initlizing containers
	l = new double[n_x - 1];
	Mu = new double[n_x - 1];
	z = new double[n_x - 1];
	b = new double[n_x - 1];	
	U = new double*[n_t];	
	for(size_t i=0;i<n_t;i++){
		U[i]=new double[n_x];			
	}	


	ImposeBoundaryCondition(n_t,n_x,U, del_t, del_x);


	//Solving for different schemes

	if (_PROBLEM == 2){
		CrankNicolson(U, n_t, n_x, 0.01, del_t, del_x,l,Mu,z,b);
	}
	else if (_PROBLEM == 3){
		CrankNicolson(U, n_t, n_x, 1.0, del_t, del_x, l, Mu, z, b);
	}
	else if (_PROBLEM == 4) {
		BDF2(U, n_t, n_x, 1.0, del_t, del_x, l, Mu, z, b);
	}

	
		
	//Plotter plot;
	//plot.ColorIsoCountoring("results.ps",del_t,n_t,del_x,n_x,U);	
	
	//Print out raw data for matlab to plot graphs
	std::string fname = "results" + std::to_string(del_t) + "_" + std::to_string(del_x) + ".txt";
	PrintOut(fname,U,n_t,n_x,del_t,del_x);

	return 0;
}

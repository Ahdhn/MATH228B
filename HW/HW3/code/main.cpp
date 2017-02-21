//Solving problem 1 in homework three - MATH 228B - Winter 2017
//By Ahmed H. Mahmoud
//The code solves the problems, and produces raw data file for matlab to plot 

#include "Plotter.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h> 
#define _USE_MATH_DEFINES
#include<math.h>
using namespace std;



Plotter _plot;//for plotting anywhere

inline void PrintOut(string outputfilename,double***temp,size_t n_x,size_t n_y,double del_x,double del_y,size_t time)
{
	//print out the temp at certain time
	//to be used for plotting in excel
	ofstream file(outputfilename,ios::out);
	file.precision(5);
	int j,i;
	double x,y;

	file<<"   ";
	for(i=0;i<n_x;i++){
		x=(i+1)*del_x;
		if(x==0.05 || x==0.1 ||x==0.15 ||x==0.2 || x==0.25 || x==0.3){
			file<<x<<" ";
		}
	}
	file<<endl;
	

	for(j=n_y-1;j>=0;j--){	
		y=(j+1)*del_y;
		file<<y<<"  ";
		for(i=0;i<n_x;i++){	
			x=(i+1)*del_x;
			if(x==0.05 || x==0.1 ||x==0.15 ||x==0.2 || x==0.25 || x==0.3){
				file<<temp[time][i][j]<<" ";
			}
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

inline void ImposeBoundaryCondition(size_t n_x, double del_x, size_t n_y, double del_y, double***U, size_t max_t)
{
	//set the initial and boundary conditions 
	//at time =0, T=To at 0< x< l && 0<y<w
	//at time>=0, T=T1 at 0<=x<=l && y=0
	//at time>=0, T=T2 at 0< y< w && x=0
	//at time>=0, T=T3 at 0<=x<=l && y=w
	//at time>=0, T=T4 at 0< y< w && x=l

	//initial conditions 
	size_t i, j, t;
	for (i = 1; i<n_x - 1; i++){
		for (j = 1; j<n_y - 1; j++){
			double x = double(i)*del_x;
			double y = double(j)*del_y;
			double xx = (x - 0.3)*(x - 0.3);
			double yy = (y - 0.4)*(y - 0.4);
			U[0][i][j] = exp(-100 * (xx + yy));
		}
	}

	for (t = 0; t<max_t; t++){
		for (i = 0; i<n_x; i++){
			U[t][i][0] = 0;
		}
	}

	for (t = 0; t<max_t; t++){
		for (j = 1; j<n_y - 1; j++){
			U[t][0][j] = 0;
		}
	}

	for (t = 0; t<max_t; t++){
		for (i = 0; i<n_x; i++){
			U[t][i][n_y - 1] = 0;
		}
	}

	for (t = 0; t<max_t; t++){
		for (j = 1; j<n_y - 1; j++){
			U[t][n_x - 1][j] = 0;
		}
	}
}

inline void Decomp (size_t n_x, double alfa, double gamma, double beta, double*l, double*Mu)
{
	//decompose the tri-dia matrix into two matrices
	//http://www.math.buffalo.edu/~pitman/courses/mth437/na2/node3.html
	//http://www.math.buffalo.edu/~pitman/courses/mth437/na2/node4.html
	size_t i;

	l[1]=alfa;

	for(i=2; i<n_x-1; i++){
		Mu[i-1]=gamma/l[i-1];
		l[i]=alfa-beta*Mu[i-1];		
	}
}
inline void Solv (size_t n_x,double*b, double beta, double*l, double*Mu, double*z, double*&X)
{
	//solve the decomposed matirces 
	
	z[1]=b[1]/l[1];
	for (size_t i=2; i<n_x-1; i++){ 
		z[i]=(b[i]-beta*z[i-1])/l[i];
	}

	X[n_x-2]=z[n_x-2];
	for (size_t i=n_x-3; i>=1;i--){
		X[i]=z[i]-Mu[i]*X[i+1];
	}	

}
inline void ADI(double***temp,size_t n_x,size_t n_y,double del_x,double del_y,double del_t,size_t max_t,size_t&steady_time)
{
	//solving using ADI
	//start with x-sweep then y-sweep
	//solving the tri-dia matrix using LU-decomposition method
	//http://www.math.buffalo.edu/~pitman/courses/mth437/na2/node2.html

	double cof_x,cof_y,alfa,beta,gamma,alfa1,beta1,gamma1,**temp_aux,*l,*Mu,*z,*b,*X,ac_error;
	size_t i,j,t,n,V;
	cof_x=(0.5*del_t)/(del_x*del_x);	
	cof_y=(0.5*del_t)/(del_y*del_y);

	alfa =1+2.0*cof_x;beta =-cof_x;gamma =-cof_x;//matrix cof's
	alfa1=1+2.0*cof_y;beta1=-cof_y;gamma1=-cof_y;//matrix cof's

	if(n_x>n_y){n=n_x;}//larger n
	else{n=n_y;}

	l =new double[n-1];//helper matrix
	Mu=new double[n-1];//helper matrix
	z =new double[n-1];//helper matrix
	b =new double[n-1];//RHS 

	X=new double[n-1];//used to copy the solution in temp_aux

	temp_aux=new double*[n_x];//holding values of intermediate time n+0.5
	for(i=0;i<n_x;i++){
		temp_aux[i]=new double[n_y];
	}
	//BC of the intermidate time step
	for(i=0;i<n_x;i++){
		temp_aux[i][0]=temp[0][i][0];
	}
	for(j=1;j<n_y-1;j++){
		temp_aux[0][j]=temp[0][0][j];
	}
	for(i=0;i<n_x;i++){
		temp_aux[i][n_y-1]=temp[0][i][n_y-1];
	}
	for(j=1;j<n_y-1;j++){
		temp_aux[n_x-1][j]=temp[0][n_x-1][j];
	}


	for(t=1;t<max_t;t++){		
		//x-sweep
		for(j=1;j<n_y-1;j++){	

			for(i=1;i<n_x-1;i++){ //setting up the matrix
				b[i]=cof_y*temp[t-1][i][j+1]+(1-2.0*cof_y)*temp[t-1][i][j]+cof_y*temp[t-1][i][j-1];
			}
			b[1]-=gamma*temp_aux[0][j];
			b[n_x-2]-=beta*temp_aux[n_x-1][j];

			Decomp(n_x,alfa,gamma,beta,l,Mu);//solving the tri-dia matrix
			Solv(n_x,b,beta,l,Mu,z,X);
			
			for(V=1;V<n_x-1;V++){//copying the solution into temp_aux
				temp_aux[V][j]=X[V];
			}
		}
		//y-sweep
		for(i=1;i<n_x-1;i++){

			for(j=1;j<n_y-1;j++){//setting up the matrix
				b[j]=cof_x*temp_aux[i+1][j]+(1-2.0*cof_x)*temp_aux[i][j]+cof_x*temp_aux[i-1][j];
			}
			b[1]-=gamma1*temp[t][i][0];
			b[n_y-2]-=beta1*temp[t][i][n_y-1];

			Decomp(n_y,alfa1,gamma1,beta1,l,Mu);//solving the tri-dia matrix
			Solv(n_y,b,beta1,l,Mu,z,temp[t][i]);
		}
		
		ac_error=0;
		for(i=1;i<n_x-1;i++){
			for(j=1;j<n_y-1;j++){
				ac_error+=abs(temp[t][i][j]-temp[t-1][i][j]);				
			}
		}
		if(ac_error<0.01){break;}//break if we reached steady state	
		if(t*del_t==40.0 || t*del_t==10.0){
			PrintOut("ADI.txt",temp,n_x,n_y,del_x,del_y,t);
			_plot.ColorIsoCountoring("ADI.ps",del_x,n_x,del_y,n_y,temp[t]);
		}
	}

	if(t==max_t){
		steady_time=0;
		cout<<"\n WARNING!! We may have not reached a steady state solution"<<endl;
		cout<<"Increase max_t (maximum number of time steps)"<<endl;
		system("pause");

	}else{
		steady_time=t-1;
		cout<<"\n Reached steady state solution after "<<steady_time<<" iterations"<<endl;
	}
	_plot.ColorIsoCountoring("ADI.ps",del_x,n_x,del_y,n_y,temp[t-1]);

}


int main (int argc, char *argv[])
{
	double***U;//solution (temperature) for n time to all grid points	
	double l, w, del_t, del_x, del_y;
	size_t n_x, n_y, max_t, i, steady_time;


	//user input values	
	cout << "Solving problem 1 in homework three - MATH 228B - Winter 2017" << endl;
	
	cout<<"\n******** GEOMETRY & BOUNDARY CONDITIONS ********"<<endl;
	cout<<"Bar length L (m) = ";//x-direction
	cin>>l;
	cout<<"Bar width W (m) = ";//y-direction
	cin>>w;	
	cout<<"Max number of cells in length direction IMAX =";
	cin>>n_x;
	cout<<"Max number of cells in width direction JMAX =";
	cin>>n_y;
	cout<<"Time step = ";
	cin>>del_t;

	del_x = l / double(n_x);
	del_y = w / double(n_y);


	max_t=1.0;
	U = new double**[2];//we only need to store two time steps 
	for (size_t t = 0; t < 2; t++){ 
		U[t] = new double*[n_x];
		for (size_t i = 0; i < n_x; i++){ 
			U[t][i] = new double[n_y];
		}
	}

	ImposeBoundaryCondition(n_x, del_x, n_y, del_y, U, max_t);
		
	

	ADI(U, n_x, n_y, del_x, del_y, del_t, max_t, steady_time);
	PrintOut("ADI.txt", U, n_x, n_y, del_x, del_y, steady_time);
	PrintOutTransAtPosition(0.1, 0.05, U, n_x, n_y, del_x, del_y, del_t, steady_time);
	PrintOutTransAtPosition(0.15, 0.1, U, n_x, n_y, del_x, del_y, del_t, steady_time);
	PrintOutTransAtPosition(0.1, 0.3, U, n_x, n_y, del_x, del_y, del_t, steady_time);



	char dummy;
	cout<<"\n Press any key to quit"<<endl;
	cin>>dummy;
	return 0;


}
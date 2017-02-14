#include "Plotter.h"
#include <iostream>
#include <string.h>  
#include <fstream>
#include <sstream>
#define _USE_MATH_DEFINES
#include <math.h>
using namespace std;

void Plotter::ColorIsoCountoring(string outputfilename,double del_x, size_t nx, double del_y, size_t ny,double**val)
{
	ofstream file(outputfilename,ios::out);	
	file << "%!PS-Adobe-3.0" << endl;
    file << "72 72 scale     % one unit = one inch" << endl;
		
	double scale_x, scale_y,shift_x, shift_y,scale,h,min_val,max_val,spc,function,u,w,tol,u1,w1,u2,w2,u3,w3,u4,w4;
	size_t num_iso,i,j,k,num,j2,k2;
    
	tol=del_x*nx*10E-9;

	scale_x = 6.5 / (del_x*(nx-1));
	scale_y = 9.0 / (del_y*(ny-1));	
	if(scale_x<scale_y){    
		scale = scale_x; 
		shift_x = 0.9;
		shift_y = (11.0-(ny-1)*del_y*scale)/2.0;
	}else{    
		scale = scale_y;    
		shift_x = (8.3-(nx-1)*del_x*scale)/2.0;
		shift_y = 1.35;
	}  
	file << shift_x << " " << shift_y << " translate" << endl;
	num_iso=80;

	h=255.0/(num_iso/4.0);

#pragma region Commands
	for(i=1;i<(num_iso/4);i++){
		file <<"/line_"<<i<<"      %stack: x1 y1_ x2 y2"<< endl;
		file << "{newpath" << endl;
		file << "moveto" << endl;
		file << "lineto" << endl;
		file << "1 "<<((i-1)*h)/255.0<<" 0 setrgbcolor" << endl;
		file << " 0.03 setlinewidth" << endl;
		file << " closepath" << endl;
		file << " stroke" << endl;
		file << "} def" << endl;
	}
	for(i=(num_iso/4);i<(num_iso/2);i++){
		file <<"/line_"<<i<<"      %stack: x1 y1_ x2 y2"<< endl;
		file << "{newpath" << endl;
		file << "moveto" << endl;
		file << "lineto" << endl;
		file << (255-(i-size_t(num_iso/4))*h)/255.0<<" 1 0 setrgbcolor" << endl;
		file << " 0.03 setlinewidth" << endl;
		file << " closepath" << endl;
		file << " stroke" << endl;
		file << "} def" << endl;
	}
	for(i=(num_iso/2);i<(3*num_iso/4);i++){
		file <<"/line_"<<i<<"      %stack: x1 y1_ x2 y2"<< endl;
		file << "{newpath" << endl;
		file << "moveto" << endl;
		file << "lineto" << endl;
		file << " 0 1 "<<((i-size_t(num_iso/2))*h)/255<<" setrgbcolor" << endl;
		file << " 0.03 setlinewidth" << endl;
		file << " closepath" << endl;
		file << " stroke" << endl;
		file << "} def" << endl;
	}	
	for(i=(3*num_iso/4);i<num_iso;i++){
		file <<"/line_"<<i<<"      %stack: x1 y1_ x2 y2"<< endl;
		file << "{newpath" << endl;
		file << "moveto" << endl;
		file << "lineto" << endl;
		file << " 0 "<< (255-(i-size_t((3*num_iso)/4))*h)/255<<" 1 setrgbcolor" << endl;
		file << " 0.03 setlinewidth" << endl;
		file << " closepath" << endl;
		file << " stroke" << endl;
		file << "} def" << endl;
	}

	file <<"/line_"<<num_iso<<"      %stack: x1 y1_ x2 y2"<< endl;
	file << "{newpath" << endl;
	file << "moveto" << endl;
	file << "lineto" << endl;
	file << " 0 0 1 setrgbcolor" << endl;
	file << " 0.03 setlinewidth" << endl;
	file << " closepath" << endl;
	file << " stroke" << endl;
    file << "} def" << endl;

	file << "/seg      % stack: x1 y1_ x2 y2" << endl;
    file << "{newpath" << endl; 
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
	file << " 0.0 0.0 0.0 setrgbcolor" << endl;
    file << " 0.007 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;

	file << "/quad_bold      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << endl;
    file << "{newpath" << endl; 
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
	file << " 0.0 0.0 0.0 setrgbcolor" << endl;
    file << " 0.002 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;

	file << "/quad      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << endl;
    file << "{newpath" << endl; 
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
	file << " fill" << endl;
    file << " 0.02 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;
#pragma endregion 

	//find max and min value

	min_val=val[0][0];
	max_val=val[0][0];

	for(i=0;i<nx;i++){
		for(j=0;j<ny;j++){
			if(max_val<val[i][j]){max_val=val[i][j];}
			if(min_val>val[i][j]){min_val=val[i][j];}		
		}
	}

	spc=(max_val-min_val)/(num_iso-1);

	for(i=0;i<num_iso;i++){
		function=max_val-i*spc;
		//cout<<"iso-contour no.["<<i<<"] val="<<val<<endl;
		
		for(j=0;j<nx-1;j++){
			u=j*del_x;
			j2=j+1;
			

			for(k=0;k<ny-1;k++){
				w=k*del_y;
				num=0;
				k2=k+1;
												
				if(function==val[j][k]&&function==val[j2][k] || function==val[j][k]&&function==val[j][k2] || function==val[j][k2]&&function==val[j2][k2] || function==val[j2][k]&&function==val[j2][k2]){					
					if(abs(function-val[j][k])<=tol && abs(function-val[j2][k])<=tol){						
						u1=u;
						w1=w;
						u2=u+del_x;
						w2=w;
						num++;
					}
					if(abs(function-val[j2][k])<=tol && abs(function-val[j2][k2])<=tol){
						if(num==0){
							u1=u+del_x;
							w1=w;
							u2=u+del_x;
							w2=w+del_y;
						}else{
							u3=u+del_x;
							w3=w+del_y;
						}
						num++;
					}
					if(abs(function-val[j][k2])<=tol && abs(function-val[j2][k2])<=tol){
						if(num==0){
							u1=u;
							w1=w+del_y;
							u2=u+del_x;
							w2=w+del_y;
						}else if(num==1){
							u3=u;
							w3=w+del_y;
						}else if(num==2){
							u4=u;
							w4=w+del_y;
						}
						num++;
					}
					if(abs(function-val[j][k])<=tol && abs(function-val[j][k2])<=tol){		
						if(num==0){
							u1=u;
							w1=w;
							u2=u;
							w2=w+del_y;					
							num++;
						}else if(num==1){
							u3=u;
							w3=w;
							u4=u;
							w4=w+del_y;					
							num+=2;
						}else if( num==2){
							cout<<"Error (0) at ColorIsoCountoring()"<<endl;							
						}						
					}					
					if(num==3){
						file<< u1*scale   <<"  "<< w1*scale <<"  ";
						file<< u2*scale   <<"  "<< w2*scale <<"  ";
						file<< u3*scale   <<"  "<< w3*scale <<"  ";
						file<< u4*scale   <<"  "<< w4*scale <<"  ";
						//file << "quad"<< endl;
						continue;

					}
					if(num>1){
						cout<<"Error(1) ColorIsoCountoring()... increase n "<<endl;						
						continue;
						//system("pause");
					}
				}else{
					if(function>val[j][k]&&function<=val[j2][k] || function<val[j][k]&&function>=val[j2][k]){	

						h=(function-val[j][k])*del_x/(val[j2][k]-val[j][k]);
						u1=u+h;
						w1=w;
						num++;
					}
					if(function>val[j][k]&&function<=val[j][k2] || function<val[j][k]&&function>=val[j][k2]){						
						h=(function-val[j][k])*del_y/(val[j][k2]-val[j][k]);
						if(num==0){
							u1=u;
							w1=w+h;
						}else{
							u2=u;
							w2=w+h;
						}
						num++;
					}
					if(function>val[j][k2]&&function<=val[j2][k2] || function<val[j][k2]&&function>=val[j2][k2]){
						
						h=(function-val[j][k2])*del_x/(val[j2][k2]-val[j][k2]);
						if(num==0){
							u1=u+h;
							w1=w+del_y;
						}else if(num==1){
							u2=u+h;
							w2=w+del_y;
						}else{
							u3=u+h;
							w3=w+del_y;
						}
						num++;
					}
					if(function>val[j2][k]&&function<=val[j2][k2] || function<val[j2][k]&&function>=val[j2][k2]){
						
						h=(function-val[j2][k])*del_y/(val[j2][k2]-val[j2][k]);
						if(num==0){
							u1=u+del_x;
							w1=w+h;
						}else if(num==1){
							u2=u+del_x;
							w2=w+h;
						}else if(num==2){
							u3=u+del_x;
							w3=w+h;
						}else{
							u4=u+del_x;
							w4=w+h;
						}
						num++;
					}
				}	

				if(num>0){
					if(num==2){
						file<< u1*scale   <<"  "<< w1*scale <<"  ";
						file<< u2*scale   <<"  "<< w2*scale <<"  ";
						file << "line_"<<i+1<< endl;
					}else if(num==3){
						file<< u1*scale   <<"  "<< w1*scale <<"  ";
						file<< u2*scale   <<"  "<< w2*scale <<"  ";
						file << "line_"<<i+1<< endl;
						file<< u2*scale   <<"  "<< w2*scale <<"  ";
						file<< u3*scale   <<"  "<< w3*scale <<"  ";
						file << "line_"<<i+1<< endl;
					}else if(num==4){
						file<< u1*scale   <<"  "<< w1*scale <<"  ";
						file<< u2*scale   <<"  "<< w2*scale <<"  ";
						file << "line_"<<i+1<< endl;
						file<< u2*scale   <<"  "<< w2*scale <<"  ";
						file<< u3*scale   <<"  "<< w3*scale <<"  ";
						file << "line_"<<i+1<< endl;
						file<< u3*scale   <<"  "<< w3*scale <<"  ";
						file<< u4*scale   <<"  "<< w4*scale <<"  ";
						file << "line_"<<i+1<< endl;
					}
				}			
			}
		}
	}


	// plot domain boundaries in bold
    file << 0.0 * scale << "  " << 0.0 * scale << "  ";  
    file << 0.0 * scale << "  " << (ny-1)*del_y * scale << "  ";          
	file << (nx-1)*del_x* scale << "  " << (ny-1)*del_y * scale << "  ";          
	file << (nx-1)*del_x* scale << "  " << 0.0 * scale << "  ";          
    file << "quad_bold"      << endl;


	//some axis libles
	file << "/Times-Roman findfont" << endl; //start of x-axis
	file << "0.15 scalefont" << endl;
	file << "setfont" << endl;
	file << -0.05*scale << " " << 0.0*scale << " moveto" << endl;
	file << "(X=0) show" << endl;

	file << "/Times-Roman findfont" << endl; //end of x-axis
	file << "0.15 scalefont" << endl;
	file << "setfont" << endl;	
	file << -0.05*scale << " " << (ny-1.2)*del_y*scale << " moveto" << endl;
	//file << "(X=L) show" << endl;
	file << "(X=" << double(ny - 1)*double(del_y) << ") show" << endl;

	file << "/Times-Roman findfont" << endl; //start of t-axis
	file << "0.15 scalefont" << endl;
	file << "setfont" << endl;
	file << 0.0*scale << " " << -0.05*scale << " moveto" << endl;
	file << "(Time=0) show" << endl;

	file << "/Times-Roman findfont" << endl; //end of t-axis	
	file << "0.15 scalefont" << endl;
	file << "setfont" << endl;
	file << (nx)*del_x*scale << " " << -0.05*scale << " moveto" << endl;
	file << "(Time="<< double(nx-1)*double(del_x)<< ") show" << endl;

	file << "showpage"<< endl;
}

	
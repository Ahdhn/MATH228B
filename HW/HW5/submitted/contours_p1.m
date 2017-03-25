%Plotting solution for problem 2 and 3 in homework one - MATH 228B - Winter 2017
%By Ahmed H. Mahmoud
%Read raw data file (fileID) and plot
clc;
%clear;
%clf;


infile_u = 'U_Step_vanLeer_.txt';
data_u = load(infile_u, '-ascii');
x_u = data_u(:,1);%col 1
y_u = data_u(:,2);%col 2


infile_anal = 'init.txt';
data_u_anal = load(infile_anal, '-ascii');
x_u_anal = data_u_anal(:,1);%col 1
y_u_anal = data_u_anal(:,2);%col 2

plot(x_u,y_u,'g');
hold on;
plot(x_u_anal,y_u_anal,'m');
xlim([0 1]);
%title('Beam-Warming');
%title('Lax-Wendroff');
%title('MC');
%title('Minmod');
%title('Superbee');
%title('Upwinding');
title('van Leer');
xlabel('X');
legend('U_{computed}', 'U_{initial}');



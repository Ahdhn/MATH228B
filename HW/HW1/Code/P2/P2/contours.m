%Plotting solution for problem 2 and 3 in homework one - MATH 228B - Winter 2017
%By Ahmed H. Mahmoud
%Read raw data file (fileID) and plot
clc;
clear;
clf;

fileID=fopen('results0.100000_0.020000.txt','r');
sz = fscanf(fileID, '%f',[1,2]);

z =fscanf(fileID, '%f',sz);

contourf(z);
set(gca,'XTickLabel',{});
set(gca,'YTickLabel',{});
title('\Delta t = 1/10, \Delta x = 1/50')
xlabel('Time  \rightarrow');
ylabel('Space \rightarrow');
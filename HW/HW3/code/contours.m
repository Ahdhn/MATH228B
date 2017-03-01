%Plotting solution for problem 2 and 3 in homework one - MATH 228B - Winter 2017
%By Ahmed H. Mahmoud
%Read raw data file (fileID) and plot
clc;
%clear;
%clf;

fileID=fopen('p2_data/part3/p2_part3_0.txt','r');
sz = fscanf(fileID, '%f',[1,2]);

z =fscanf(fileID, '%f',sz)';

figID = contourf(z,50);
zlim([-3 3]);

set(gca, 'YDir', 'reverse');
shading interp;
set(gca,'XTickLabel',{});
set(gca,'YTickLabel',{});
colorbar();
%title('time = 0.39 sec ')
fclose(fileID);

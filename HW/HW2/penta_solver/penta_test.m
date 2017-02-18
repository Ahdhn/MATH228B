% This scipt was written to test penta solver for
% MATH 501 computer assignment  CSUF
%by Nasser M. Abbasi. Feb 27,2007

clear all;
clc;
format long;

test=1;
fprintf('==========>test %d\n',test);
% A=[1 8 1 0 0 0 0 0 0;
%    1 8 4 1 0 0 0 0 0;   
%    1 8 1 1 1 0 0 0 0;   
%    0 1 1 6 4 1 0 0 0;   
%    0 0 1 6 4 5 1 0 0;   
%    0 0 0 1 3 5 1 1 0;   
%    0 0 0 0 1 1 4 8 1;   
%    0 0 0 0 0 1 4 8 9;   
%    0 0 0 0 0 0 1 9 4]
%A=[1 8 0 6 0 0 0 0 0;
%   1 8 4 0 4 0 0 0 0;   
%   0 8 1 0 0 1 0 0 0;   
%   1 0 0 6 4 0 5 0 0;   
%   0 8 0 6 4 5 0 5 0;   
%   0 0 4 0 3 5 0 0 1;   
%   0 0 0 6 0 0 4 8 0;   
%   0 0 0 0 1 0 4 8 9;   
%   0 0 0 0 0 4 0 9 4]
%A=[15   8 -6   0   0   0;
%   -2  12 -4  -4   0   0;
%   -6  -4 19  -9   4   0;
%    0  -1 -9  21   6   7; 
%    0   0  9  10  11   8;
%    0   0  0  10  -2   2]
A = zeros(25,25);
for i=1:1:25    
    %diagonal
    A(i,i)=1.0002;
end

%first sub and super diagonal 
for i=1:1:25
    for j=1:1:25
        
        if (i-j)==1 && mod(j,5)>0
            A(i,j) = -0.00005;
        end
        
        if(j-i)==1 && mod(j-1,5)>0
            A(i,j) = -0.00005;
        end
    end    
end
%outer sub and super diagonal 
for i=1:1:25
    if(i+5)<25
        A(i,i+5) = -0.00005;
    end
     if(i-5)>0
        A(i,i-5) = -0.00005;
    end
    
end
b =[0.000735644 0.108378183 0.062185974 0.000141031 8.08E-09 0.003893875 0.573661329 0.329159319 0.000746499 4.28E-08 7.99E-05 0.011770308 0.006753647 1.53E-05 8.76E-10 1.03E-08 1.52E-06 8.70E-07 1.96E-09 7.57E-14 3.15E-13 4.67E-11 2.68E-11 5.95E-14 5.28E-19];

fprintf('=======>Matlab result\n');
A\b'
fprintf('=======>our result\n');
x=nma_pentaSolve(A,b)


test=test+1;
fprintf('==========>test %d\n',test);
A=[15  -2 -6    0;
   -2  12 -4   -4;
   -6  -4 19   -9;
   0   -1 -9   21]
b=[300 0 0 0]
fprintf('=======>our result\n');
x=nma_pentaSolve(A,b)
fprintf('=======>Matlab result\n');
A\b'


test=test+1;
fprintf('==========>test %d\n',test);
A=[15   8 -6   0   0   0  0;
   -2  12 -4  -4   0   0  0;
   -6  -4 19  -9   4   0  0;
    0  -1 -9  21   6   7  0; 
    0   0  9  10  11   8  3;
    0   0  0  10  -2   2  4;
    0   0  0  0   -2   2  4]
b=[300 0 0 0 1 2 6]
fprintf('=======>our result\n');
x=nma_pentaSolve(A,b)
fprintf('=======>Matlab result\n');
A\b'

test=test+1;
fprintf('==========>test %d\n',test);
A=[15   8;
   0  12]
b=[300 0]
fprintf('=======>our result\n');
x=nma_pentaSolve(A,b)
fprintf('=======>Matlab result\n');
A\b'


test=test+1;
fprintf('==========>test %d\n',test);
A=[15   8 -6   0   0   0  0  0;
   -2  12 -4  -4   0   0  0  0;
   -6  -4 19  -9   4   0  0  0;
    0  -1 -9  21   6   7  0  0; 
    0   0  9  10  11   8  3  0;
    0   0  0  10  -2   2  4  3;
    0   0  0  0   -2   2  4  7;
    0   0  0  0   0    4  8  9]
b=[300 0 0 0 1 2 6 10]
fprintf('=======>our result\n');
x=nma_pentaSolve(A,b)
fprintf('=======>Matlab result\n');
A\b'



function [A,b]=nma_gaussian_elimination(A,b)
%function [U,b_new]=nma_gaussian_elimination(A,b)
%
%Perform Gaussian Elimination on a matrix and the
%updated the corresponding b matrix as well.
%
%This function is used as part of the routines I written
%during work on Math 501, at CSFU, winter 2007
%
%INPUT
%  A: matrix. does not have to be square
%  b: vector. from Ax=b
%
%OUTPUT:
%  U: upper triangular matrix
%  b_new: updated b vector
%
%by Nasser Abbasi. Feb 27,2007

    if nargin ~=2
        error('2 arguments required');
    end

    if ~isnumeric(A)|~isnumeric(b)
        error('A and b must be numeric');
    end

    [nRow,nCol]=size(A);
    [b_nRow,b_nCol]=size(b);
    if b_nCol>1
        error('b must be a vector');
    end

    if b_nRow~=nRow
        error('b must be same size as A');
    end

    pivot=A(1,1);
    [nRow,nCol]=size(A);
    for i=2:nRow
        multiplier=A(i,1)/pivot;
        b(i)=b(i)-multiplier*b(1);
        for j=1:nCol
            A(i,j)=A(i,j)-multiplier*A(1,j);
        end
    end
   
end
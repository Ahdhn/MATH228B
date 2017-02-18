function x=nma_pentaSolve(A,b)
%function x=nma_pentaSolve(A,b)
%Solves a penta-diagonal Ax=b system
%
%INPUT:
%  A  an nxn Matrix
%  b  an vector of length n
%
%OUTPUT
%  x  vector of length n, the solution for Ax=b
%
%Algorithm Overview:
% The matrix is banded matrix. Due to this, saving in processing
% is achived by only processing the non-zero elements along the band.
% see full report for more details
%by Nasser M. Abbasi. Feb 27,2007

%March 7, 2007

    if nargin ~=2
        error('2 arguments required');
    end

    if ~isnumeric(A)|~isnumeric(b)
        error('A and b must be numeric');
    end

    [nRow,nCol]=size(A);
    if nRow~=nCol
        error('A must be square');
    end

    b=b(:);

    [b_nRow,b_nCol]=size(b);
    if b_nCol>1
        error('b must be a vector');
    end

    if b_nRow~=nRow
        error('b must be same size as A');
    end

    [U,b_new]=penta_elimination(A,b);
    x=nma_pentaBackSub(U,b_new);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The algorithm checks for the band width, and process elements
% within the band as using standard G.E. This saves processing 
% time compared with the G.E. which process the whole matrix.
%
% Since the matrix is sparse, we save time by processing only the
% non-zero band.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,b]=penta_elimination(A,b)

    TRUE=1;
    FALSE=0;

    [nRow,nCol]=size(A);
    pivot=1;
    more_bands=TRUE;

    while more_bands
        r=find(abs(A(pivot:end,pivot))>4*eps);
        r=r+(pivot-1);  %find is relative to current band, so adjust to A
        end_row=r(end);
        r=find(abs(A(end,:))>4*eps);
        end_column=r(end);

        [U,b_new]=nma_gaussian_elimination...
            (A(pivot:end_row,pivot:end_column),b(pivot:end_row));
        A(pivot:end_row,pivot:end_column)=U;
        b(pivot:end_row)=b_new;

        if pivot==nRow
            more_bands=FALSE;
        else
            pivot=pivot+1;
        end
    end
    
end


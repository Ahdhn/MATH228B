function x=nma_pentaBackSub(U,b)
%function x=nma_pentaBackSub(U,b)
%
%Does backsub on an upper diagonal U type matrix
%special backsub in that this function can
%detected banded U and will only process the non-zero
%entries in its backsub process
%
%INPUT:
%  U: special upper triangular matrix. Can be the
%     result of doing special Gaussian elimination on penta-diagonal
%  b: the rhs in the Ux=b
%
%OUTPUT: the solution to Ux=b
%
%by Nasser M. Abbasi. Feb 27,2007

    if nargin ~=2
        error('2 arguments required');
    end

    if ~isnumeric(U)|~isnumeric(b)
        error('U and b must be numeric');
    end

    [nRow,nCol]=size(U);
    [b_nRow,b_nCol]=size(b);
    if b_nCol>1
        error('b must be a vector');
    end

    if b_nRow~=nRow
        error('b must be same size as U');
    end

    [nRow,nCol]=size(U);
    n=nRow;
    x=zeros(n,1);
    x(n)=b(n)/U(n,n);

    for i=n-1:-1:1
        r=find(abs(U(i,:))>4*eps);
        row=U(i,r(1):r(end));
        x(i)=(b(i)-row(2:end)*x(i+1:i+length(row)-1))/row(1);
    end

end
%by Nasser Abbasi. Feb 27,2007
function x=nma_backsolve(A,b)

[n,nCols]=size(A);

x=zeros(n,1);

x(n)=b(n)/A(n,n);
x(n-1)=(b(n-1)-A(n-1,n)*x(n))/A(n-1,n-1); 

for k=n-2:-1:1
   x(n-2)=(   b(n-2)-  A(n-2,n-1)*x(n-1) - A(n-2,n)*x(n)  )/A(n-2,n-2);   
end


end
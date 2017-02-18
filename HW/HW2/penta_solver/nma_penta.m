%by Nasser M. Abbasi. Feb 27,2007
function [A,b]=nma_penta(A,b)
[nRows,nCols]=size(A);

d=diag(A);

s=zeros(nRows-1,1);
for i=2:nRows
    s(i-1)=A(i,i-1);
end

u=zeros(nRows-2,1);
for i=3:nRows
    u(i-2)=A(i,i-2);
end

t=zeros(nRows-1,1);
for i=1:nRows-1
    t(i)=A(i,i+1);
end

v=zeros(nRows-2,1);
for i=1:nRows-2
    v(i)=A(i,i+2);
end


j=0;
for k=2:nRows
   j=j+1;
   m=s(j)/d(j);
   b(k)=b(k)-m*b(k-1);
   
   A(k,k-1)=A(k,k-1)-m*A(k-1,k-1);
   
   s(j)=A(k,k-1);
   
   A(k,k)=d(j+1)-m*t(j);
   d(j+1)=A(k,k);
   
   if k+1<=nRows
     A(k,k+1)=t(j+1)-m*v(j);
     t(j+1)=A(k,k+1);
   end
   
   if k+1<=nRows
     m=u(j)/d(j);
     A(k+1,k-1)=A(k+1,k-1)-m*A(k,k-1);
     u(j)=A(k+1,k-1);
     
     A(k+1,k)=s(j+1)-m*d(j+1);
     s(j+1)=A(k+1,k);
     
     A(k+1,k+1)=d(j+2)-m*t(j+1);
     d(j+2)=A(k+1,k+1);
     
     if k+2<nRows
       A(k+1,k+2)=t(j+2)-m*v(j+1);
       t(j+2)=A(k+1,k+2);
     end
     
   end

end


end
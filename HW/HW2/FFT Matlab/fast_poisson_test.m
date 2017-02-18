%
% test the 2D Poisson solver by comparing with a direct solve
%   of the problem (a-b*L)u=f using random a,b, and f
%
fprintf('%6s%12s%12s%12s%12s \n','nx','time (ge)','time (fft)','ratio','max diff');
for q = 5:11
  N = 2^q-1;
  L = 1;
  dx = L/(N+1);

  % pick random parameters and random right-hand size
  %
  a = rand;
  b = rand;
  b = b/dx^2;
  f = rand(N,N);
  
  % form the matrix to invert
  %
  L = lap2d(N,N);
  A = a*speye(N*N) - b*L;
  
  tic;
  u1 = A\f(:);
  t1 = toc;
  
  tic;
  u2 =  fast_poisson_2D(f,a,b);
  t2=toc;
  
  fprintf('%6i%12.2e%12.2e%12f%12.2e \n',...
      N,t1,t2,t1/t2,max(abs(u1(:)-u2(:))));
  
end  
  






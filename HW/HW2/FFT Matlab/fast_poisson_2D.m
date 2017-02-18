%
% solve the equation 
%   (a*I-b*L)u=f 
% where L corresponds to a 2D discrete Laplacian with dx=1 on an equal
% spaced node-centered grid, assuming Dirichlet boundary conditions
% 
%  f is an nx by ny array of the right hand side
%  a,b are scalars
%
%
function u = fast_poisson_2D(f,a,b);
  % record the size
  %
  [nx,ny]=size(f);
  
  % make a matix of eigenvalues
  %
  [kx,ky]=ndgrid(1:nx,1:ny);
  L = 2*(cos(kx*pi/(nx+1))-1) + 2*(cos(ky*pi/(ny+1))-1);
  L = a - b*L;

  % perform the solve
  %
  fhat = dst2( f );
  uhat = fhat./L;
  u    = idst2( uhat );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
function fhat =dst2(f)
  f1 = dst(f);
  fhat = dst(f1.');
  fhat = fhat.';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

function f = idst2(fhat);
  f1 = idst( fhat );
  f  = idst( f1.' );
  f  = f.';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
function fhat = dst(f);
  
  % record the size
  %
  [nx,ny]=size(f);
  
  % expand f to be an odd function
  %
  fext = [zeros(1,ny); f; zeros(1,ny); -flipud(f)];
  
  % fourier transform the x-direction
  %
  fexthat = fft(fext);
  
  % extract the sine transform
  %
  fhat = imag( fexthat(2:(nx+1),:) );
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = idst(fhat);
  
  % record the size
  %
  [nx,ny]=size(fhat);

  % expand fhat
  %
  i=sqrt(-1);
  fexthat = [zeros(1,ny); i*fhat; zeros(1,ny);  -i*flipud(fhat)];
  
  % back to real space
  %
  fext = ifft( fexthat );
  f    = real( fext((2:(nx+1)),:) );
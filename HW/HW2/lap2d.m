%
% lap2d.m
%
%   form the (scaled) matrix for the 2D Laplacian for Dirichlet boundary
%   conditions on a rectangular node-centered nx by ny grid
%   
%   input:  n -- number of grid points in x-direction (no bdy pts)
%           m -- number of grid points in y-direction
%
%   output: L2 -- (n*m) x (n*m) sparse matrix for discrete Laplacian
%
function L2 = lap2d(nx,ny);
    
    % make 1D Laplacians
    %
    Lx = lap1d(nx);
    Ly = lap1d(ny);
    
    % make 1D identities
    %
    Ix = speye(nx);
    Iy = speye(ny);
    
    % form 2D matrix from kron
    %
    L2 = kron(Iy,Lx) + kron(Ly,Ix);
        
%
% function: lap1d -- form the (scaled) 1D Laplacian for Dirichlet
%                    boundary conditions on a node-centered grid
%
%   input:  n -- number of grid points (no bdy pts)
%
%   output: L -- n x n sparse matrix for discrete Laplacian
%
function L = lap1d(n)
    e = ones(n,1);
    L = spdiags([ e -2*e e], [-1 0 1],n,n);
    

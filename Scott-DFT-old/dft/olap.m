function [S] =  olap(psi)

% [N_g, N_o] = size(psi);
n = size(psi);
N_g = prod(n(1:3));
psi_r = reshape(psi,N_g,n(4));
S = psi_r' * psi_r;

% Orthogonalize the plane wave basis
%  Compute Overlap matrix S    [ N_o^2 * N_g , (N_o x N_o) matrix]
%       S = PSI^t * PSI 
%   diagonalize S -> Lambda   [ N_o^3, QR repeatedly)
%   Apply Lambda to PSI       
%
%   PSI is N_g x N_o

function [psiO] =  ortho(psi)

% [N_g, N_o] = size(psi);
n = size(psi);
N_g = prod(n(1:3));
psi_u = reshape(psi,N_g,n(4));
psiO = cgrscho(psi_u);
%psiO = qr(psi_u);
psiO = reshape(psiO,n);

% S = psi_u' * psi_u;
% [A,e] = eig(S);
% psiO = A * psi_u;



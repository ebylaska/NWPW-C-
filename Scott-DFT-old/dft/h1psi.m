%   PSI is N_g x N_o (Complex)
%   N_g is m x m x m

function [h1psi_G] =  h1psi(psi,Vloc_G,G,omega)

C =  -(3.0/(4*pi^2))^(1.0/3.0);

% Needed because matlab scales the inverse FFT
s=size(psi); scale = prod(s(1:3));

psi_r = zeros(size(psi));
N_o = s(4);
for k = 1 : N_o
   psi_r(:,:,:,k) = scale*ifftn(squeeze(psi(:,:,:,k)));
end
   rho_r = 2.0*sum(conj(psi_r) .* psi_r,4)/omega;
   if (~isreal(rho_r))
      'Is not real!' 
   end
   %display([ 'N: '  num2str( sum(sum(sum(rho_r)))*omega/scale) ]);
   Vxc_r = C * rho_r .^ (1.0/3.0);

   rho_G = fftn(rho_r)/scale;
   G2 = squeeze(dot(G,G,1));
% Project against zero divide
   G2(1,1,1) = 1.0;
   VH_G = (4*pi)/(G2) .* rho_G;
   VH_G(1,1,1) = 0.0;

   %display([ 'E_H: '  num2str(omega*0.5*sum(sum(sum(conj(rho_G) .* VH_G))))])
   VH_r = scale*ifftn(VH_G+Vloc_G/omega);
for k = 1 : N_o
% Vloc_r is in VH_r
   h1psi_r = (VH_r + Vxc_r) .* squeeze(psi_r(:,:,:,k));
   h1psi_G(:,:,:,k) = fftn(h1psi_r)/scale;
end

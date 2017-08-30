s=size(psi); scale = prod(s(1:3));

psi_r = zeros(size(psi));
N_o = s(4);
for k = 1 : N_o
   psi_r(:,:,:,k) = scale*ifftn(squeeze(psi(:,:,:,k)));
end

%   PSI is N_g x N_o
%   N_g is m x m x m

function [H2_psi] =  h2psi(psi,G)

  g2 = squeeze(dot(G,G,1));
  H2_psi = zeros(size(psi));
  for i = 1 : size(psi,4)
    H2_psi(:,:,:,i)  = (g2 .* squeeze(psi(:,:,:,i ))) * 0.5;
  end

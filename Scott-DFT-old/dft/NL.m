%%  1.0 The Non-local pseudopotential
%%  --------------------------------
%%  
%%  90% to 98% of the overall
%%  Very expensive, but very parallel
%%  Goes in parallel with the other work
%%  N_g * (N_proj * N_A)  -->  Nproj * Na ~= N_o
%%  
%%  Compute for each orbital i
%%      V_NL psi_i = SUM_(I,p) | PHI_p^I> C_p <PHI_p^I PSI_i >
%%  
%%  PSI_i  : N_g x 1 column vector
%%  
%%  
%%  
%%  for  i = 0 : N_o - 1
%%     for I = 0 : N_A - 1
%%          for p = 0 : N_proj - 1
%%             Sum = 0
%%             for G = 0 : N_g -1
%%                sum += PSI_i(G) * PHI_p^I(G)
%%             Sum = Sum * C_p
%%  
%%             for G = 0 : N_g -1
%%                PSI_i(G) += PHI_p^I(G) * Sum
%%          end for
%%      end for
%%  end for
%%  
%%  
%%  For each i :
%%      Column vector times C_p * Dot product of a row vector and column vector
%%  
%%  Or, generate an N_g x N_o matrix of PSIs
%%  
%%      (N_g x 1 column vector for PHI)  *  C_p *
%%          ( 1 x N_g  row vector for PHI * N_g x N_o PSI)
%%  
%%  
%%  This can be structured into big dgemms
%%  

function [V_nl] =  NL(psi)

psi_r = psi;
for i = 1 : N_o
   psi_r(:,:,:,i) = fft(psi(:,:,:,i));
end
   rho_r = sum(psi_r .* psi_r,4);
   Vxc_r = C * rho_r^(4/3);
   rho_G = fft(rho_r);
   VH_G = (4*pi)/(G^2)* rho_G;
   VH_r = fft(VH_G);
for i = 1 : N_o
   Hpsi_r  = (VH_r .+ Vxc_r) .* psi_r(:,:,:,i);
   Hpsi_G(:,:,:,i) = fft(Hpsi_r(:,:,:,i));
end

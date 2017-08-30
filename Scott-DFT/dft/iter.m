function [psi] =  iter(alpha,L,N,N_o,N_a,N_it)
G = GGen(L,N);
psi = rand(N,N,N,N_o);
psi=ortho(psi);

n = size(psi); N_g = prod(n(1:3));

% Atomic positions - eventually we'll use realistic  data
R = rand(3,N_a);
N_a = 1;
R = zeros(3,1);

% Constants for the local potential of Silicon centered at the origin
zv =  4;
rLoc =  0.44;
C1 =  -7.336103; C2 =  0.0; C3 =  0.0; C4 = 0.0;

Vloc0 = genVloc(G, zv, rLoc, C1, C2, C3, C4);


Vloc = zeros(size(Vloc0));
for k = 1 : N_a
Vloc = Vloc +  StructFact(G,R(:,k)) .* Vloc0;
end

A = L * eye(3);
a1 = A(:,1); a2 = A(:,2); a3 = A(:,3);
omega =  dot(cross(a1,a2),a3);
% Vloc = Vloc /  omega;

% ------- end of gLocal


et = zeros(1,N_it);
et1 = zeros(1,N_it);
et2 = zeros(1,N_it);
% E_exct = 2 * (pi/L)^2
for it = 1:N_it
  h2_psiG = h2psi(psi,G);
  h1_psiG =  h1psi(psi,Vloc,G,omega);
  H_psi = h1_psiG + h2_psiG;
  E = sum(sum(sum(sum(real(conj(psi) .* H_psi)))));
  E1 = (sum(sum(sum(real(conj(squeeze(psi(:,:,:,1))) .* squeeze(H_psi(:,:,:,1)))))));
  E2 = (sum(sum(sum(real(conj(squeeze(psi(:,:,:,2))) .* squeeze(H_psi(:,:,:,2)))))));
  psi = psi - alpha * H_psi;
  psi=ortho(psi);
%  olap(psi)
 if (mod(it,10) == 0)
%    disp([ 'E: ' num2str(E) ' [' num2str(E_exct) ']' ]);
    disp([ 'E1, E2: ' num2str(E1) ', ' num2str(E2) ]);
 end
  et(it) = E;
  et1(it) = E1;
  et2(it) = E2;
%  if (abs(E - E_exct) < 1.0e-6)
%      et=et(1:i);
%      disp(['Yay after ' num2str(i) ' iterations! ']);
%      break;
%   end
end
plot(1:length(et),et,1:length(et),et1,1:length(et),et2)

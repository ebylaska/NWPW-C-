%   PSI is N_g x N_o
%   N_g is m x m x m

function [G] =  GGen(L,N)

% Simplifying assumption that the box  is a cube
% A can correspond to any box
A = L * eye(3);
a1 = A(:,1); a2 = A(:,2); a3 = A(:,3);
omega =  dot(cross(a1,a2),a3);
tpo = 2*pi/omega;

b1 = tpo * cross(a2,a3);
b2 = tpo * cross(a3,a1);
b3 = tpo * cross(a1,a2);
[I,J,K] = meshgrid(-N/2+1:N/2, -N/2+1:N/2, -N/2+1:N/2);
G = zeros(3,N,N,N);
G(1,:,:,:) = I*b1(1) + J*b2(1) + K*b3(1);
G(2,:,:,:) = I*b1(2) + J*b2(2) + K*b3(2);
G(3,:,:,:) = I*b1(3) + J*b2(3) + K*b3(3);
svec = [0 N/2+1 N/2+1 N/2+1];
G = circshift(G,svec);

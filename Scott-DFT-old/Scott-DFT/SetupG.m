% [G] = SetupG(N,a)
%
% Computes the Fourier Vectors 
%
%  Entry - 
%	N = [n1 n2 n3] - size of grid in 3 dimensions
% 	a(3,3) - lattice vectors
%
%  Exit -
%	G(3,n1,n2,n2) - Fourier vectors
%

function [G] = SetupG(N,a)
   a1 = a(:,1);
   a2 = a(:,2);
   a3 = a(:,3);
   omega = dot(cross(a1,a2),a3);
   tpo = 2*pi/omega;
   b1 =  tpo*cross(a2,a3);
   b2 =  tpo*cross(a3,a1);
   b3 =  tpo*cross(a1,a2);

   nx = N(1);
   ny = N(2);
   nz = N(3);
   nxh = nx/2;
   nyh = ny/2;
   nzh = nz/2;

   G = zeros(3,nx,ny,nz);
   for k=-nzh+1:nzh
   for j=-nyh+1:nyh
   for i=-nxh+1:nxh
     ii = i;
     jj = j;
     kk = k;
     if (ii<0) 
        ii = ii + nx;
     endif
     if (jj<0) 
        jj = jj + ny;
     endif
     if (kk<0) 
        kk = kk + nz;
     endif
     G(:,ii+1,jj+1,kk+1) = i*b1 + j*b2 + k*b3;

   endfor
   endfor
   endfor

endfunction


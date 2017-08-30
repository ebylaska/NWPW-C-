% Compute V_local
function [vloc] =  genVloc(G, zv, rLoc, C1, C2, C3, C4)
% constants: zv, rLoc, c1, c2, c3, c4
Q2 = squeeze(dot(G,G,1));
Q = sqrt(Q2);
b0 = -4*pi*zv;
b1 = sqrt((8*pi^3)*rLoc^3);
x  = Q*rLoc;
x2 = x .* x;
x4 = x2 .* x2;
x6 = x4 .* x2;
a  = exp(x2/2);
% Project against zero divide
Q2(1,1,1) = 1.0;
vloc = b0 * (a ./ Q2) + (b1 * a ) .* ( C1 ...
                                     + C2 * (3-x2)  ...
                                     + C3 * (15-10*x2 + x6)  ...
                                     + C4 * (105 - 105*x2 + 21*x4 - x6));
% Correct location (1,1,1)
vloc(1,1,1) = b1  * (C1 + 3*C2 + 15*C3 + 105*C4 ) - b0 *0.5* rLoc^2;

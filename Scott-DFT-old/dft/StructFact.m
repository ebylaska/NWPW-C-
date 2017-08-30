function [sfact] =  StructFact(G,r3)
gr = squeeze(G(1,:,:,:)) * r3(1) + squeeze(G(2,:,:,:)) * r3(2) + squeeze(G(2,:,:,:)) * r3(3);
gr = squeeze(G(1,:,:,:) * r3(1) + G(2,:,:,:) * r3(2) + G(2,:,:,:) * r3(3));
sfact = exp(-i * gr );   % Shape is shape(G)(2:4)

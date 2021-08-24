function angles = randDirAngleMJ(n, r, nsim)
% randDirAngleMJ   calculate random direction angle
%   Simulate the angles between a random direction and a fixed rank r 
%   subspace in R^n (in this case the subspace of the first r dimensions).
%
% Inputs:
%   n - dimension of vector space
%   r - dimension of subspace
%   nsim - number of simulation samples
%
% Outputs:
%   angles - nsim x 1 simulated random direction angles
%
%   Copyright (c)  Meilei Jiang 2018

    angles = zeros(nsim, 1);
    for i = 1:nsim
        vec = randn(n, 1);
        angles(i) = acosd(norm(vec(1:r))/norm(vec));
    end
end


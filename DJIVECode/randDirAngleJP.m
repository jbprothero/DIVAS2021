function angles = randDirAngleJP(n, r, nsim, rRand)
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

if nargin < 4
    rRand = 1 ;
end

    angles = zeros(nsim, rRand);
    for i = 1:nsim
        vec = randn(n, rRand);
        angles(i,:) = acosd(svd([eye(r) zeros(r, n-r)] * orth(vec)));
    end
end


function maxAngles = maxPerturbAngleMJ(validPC, VVHatCache)
% maxPerturbAngleMJ   Caluclate maximum perturbation angles
%   Detailed explanation goes here
%
% Inputs:
%   validPC - vector of logic of valid Principal Components
%   VVHatCache - array of bootstrap samples of V' * VHat
%
% Outputs:
%   maxAngles - vector of bootstrap sample of max perturbation angles
%
%   Copyright (c)  Meilei Jiang 2018

    nsim = length(VVHatCache);
    maxAngles = 90 * ones(nsim, 1);
    parfor i = 1:nsim
        VVBar = diag(validPC) * VVHatCache{i} * diag(validPC);
        maxAngles(i) = acosd(min(abs(svds(VVBar, 1, 'smallestnz')), 1));
    end
end


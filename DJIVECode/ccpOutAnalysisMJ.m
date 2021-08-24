function angleHats = ccpOutAnalysisMJ(cache_v, VBars)
% ccpOutAnalysisMJ   Calculate projected angle of each opt-vs
%   Detailed explanation goes here
%
% Inputs:
%   cache_v - cached opt-vs in the each iteration of optmization
%   VBars - cell array of adjusted signal row spaces
%
% Outputs:
%   angleHats - estimated projection angle of opt-v on each data blocks.
%
%   Copyright (c)  Meilei Jiang 2018

    nb = length(VBars);
    T = length(cache_v);
    angleHats = cell(nb, 1);
    for ib = 1:nb
        angles = -1*ones(1, T);
        for t = 1:T
            angles(t) = projAngleMJ(cache_v{t}, VBars{ib});
        end
        angleHats{ib} = angles;
    end
end


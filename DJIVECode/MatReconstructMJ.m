function [matBlockMap, matLoadingMap] = MatReconstructMJ(X, matJointV, ...
    matJointOrder, matJointRanks)
% MatReconstructMJ   Summary of this function goes here
%   Detailed explanation goes here
%
% Inputs:
%   X - d x n data matrix
%   matJointV - n x sum(r_t) joint structure basis matrix
%   matJointOrder - n x 1 cell array which keeps the order of joint blocks
%                   in matJointV
%   matJointRanks - n x 1 array which contains the rank of each joint 
%                   blocks in matJointV.
%
% Outputs:
%   matBlockMap - map between joint block index and joint block mat in X
%   matLoadingMap - map between joint block index and joint block loadings
%                   in X
%
%   Copyright (c)  Meilei Jiang 2018

    matLoading = linsolve(matJointV, X')'; %This is where the magic happens
    
    matBlockMap = containers.Map();
    matLoadingMap = containers.Map();
    cum = 0;
    for j = 1:length(matJointOrder)
        t = matJointOrder{j};
        r = matJointRanks(j);
        Vhat = matJointV(:, (cum+1):(cum+r));
        Bhat = matLoading(:, (cum+1):(cum+r));
        matBlockMap(t) = Bhat * Vhat';
        matLoadingMap(t) = Bhat;
        cum = cum + r;
    end
end


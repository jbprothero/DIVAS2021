function [VBars, phiBars, EHats, rBars]= DJIVESignalExtractMJ(datablock, ...
    dataname, nsim)
% DJIVESignalExtractMJ  Singmal matrix extraction
%      First Step of DJIVE algorithm
%
% Inputs:
%   datablock - cell array of d_k x n data matrix
%   dataname - cell array of string of each data matrix's name
%   nsim - number of bootstrap samples
%
% Outputs:
%   VBars - cell array of adjusted signal row spaces
%   phiBars - vector of perturbation angle for each data matrix
%   EHats - cell array of estimated noise matrices
%   rBars - vector of adjusted signal ranks for each data matrix
%
%   Copyright (c)  Meilei Jiang 2018
    nb = length(datablock);
    VBars = cell(nb, 1);
    EHats = cell(nb, 1);
    phiBars = 90 * ones(nb, 1);
    rBars = zeros(nb, 1);
    for ib = 1:nb
        fprintf('Signal estimation for %s\n', dataname{ib});
        [VBars{ib}, phiBars(ib), rBars(ib), EHats{ib}] = ...
            MatSignalExtractMJ(datablock{ib}, dataname{ib}, nsim);
    end        
end
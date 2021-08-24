function [outMap, keyIdxMap, anglesMap, jointBlockOrder] = DJIVEJointStrucEstimateJP(...
    VBars, UBars, phiBars, psiBars, rBars, dataname, theta0, optArgin, iprint, figdir)
% DJIVEJointStrucEstimateMJ   Estimate Joint and Paritally Joint structure
%   Estabilish a DC programming probelm to estimate each partially joint
%   structure. Using Penalty CCP algorithm to solve the DC programming 
%   problem.
%
% Inputs:
%   VBars - cell array of adjusted signal row spaces
%   phiBars - vector of perturbation angle for each data matrix
%   rBars - vector of adjusted signal ranks for each data matrix
%   dataname - a nb x 1 cell array of string of data matrix's name
%
% Outputs:
%   outMap - mapping between joint block index set and estimated partially 
%            shared structure.
%   keyIdxMap - mapping between joint block index and data blocks
%   jointBlockOrder - a string cell array to record the 
%   
%
%   Copyright (c)  Meilei Jiang 2018

    if ~exist('theta0', 'var')
        theta0 = 45;
    end
    if ~exist('optArgin', 'var')
        optArgin = [];
    end
    
    if ~exist('iprint', 'var')
        iprint = 0;
    end
    
    if ~exist('figdir', 'var')
        figdir = '';
    end

    nb = length(VBars);
    allIdx = 1:1:nb;
    curRanks = zeros(nb, 1);
    outMap = containers.Map();
    keyIdxMap = containers.Map();
    anglesMap = containers.Map();
    jointBlockOrder = {};
    flag = false;
    
    for len = nb:-1:1
        lenIdces = nchoosek(allIdx, len);
        nlen = size(lenIdces, 1);
        for i = 1:nlen
            blockIdx = lenIdces(i, :);
            blockIn = ismember(allIdx, blockIdx);
            if any(curRanks + blockIn' > rBars)
                continue;
            end
            [Vi, curRanks, angles] = BlockJointStrucEstimateJP(blockIn, dataname, ...
                VBars, phiBars, rBars, curRanks, outMap, theta0, optArgin, ...
                iprint, figdir);
            if size(Vi, 2) > 0
                t = Idx2numMJ(blockIn);
                outMap(num2str(t)) = Vi;
                keyIdxMap(num2str(t)) = blockIdx;
                anglesMap(num2str(t)) = angles;
                jointBlockOrder{end + 1} = num2str(t);
            end
            if curRanks == rBars
                fprintf('There is no room for next joint block. Stop seraching.\n')
                flag = true;
                break
            end
        end
        
        if flag
            break
        end
    end
    
end


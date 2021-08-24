function outstruct = DJIVEReconstructMJ(datablock, dataname, outMap, ...
        keyIdxMap, jointBlockOrder, doubleCenter)
% DJIVEReconstructMJ   Summary of this function goes here
%   Detailed explanation goes here
%
% Inputs:
%   datablock - cell array of d_k x n data matrices
%   dataname  - a nb x 1 cell array of string of data matrix's name
%   outMap - mapping between block index set and estimated partially shared
%            structure.
%   keyIdxMap - mapping between joint block index and data blocks
%   jointBlockOrder - a string cell array to record the 
%
% Outputs:
%   outstruct
%
%   Copyright (c)  Meilei Jiang 2018

    nb = length(datablock);
    matJointV = cell(nb, 1);
    matJointOrder = cell(nb, 1);
    matJointRanks = cell(nb, 1);
    for ib = 1:nb
        matJointV{ib} = [];
        matJointOrder{ib} = {};
        matJointRanks{ib} = [];
    end
    
    % collect joint block basis for each data matrix
    for i = 1:length(jointBlockOrder)
        t = jointBlockOrder{i};
        V = outMap(t);
        r = size(V, 2);
        blockIdx = keyIdxMap(t);
        fprintf(strcat('The rank of joint block among ', ...
            strjoin(dataname(blockIdx), ', '), ...
            ': %d \n'), r)
        for ib = blockIdx
            matJointV{ib} = [matJointV{ib}, V];
            matJointOrder{ib}{end + 1} = t;
            matJointRanks{ib}(end + 1) = r;
        end
    end
    % calculate the loadings and the joint block in each data matrix
    matLoadings = cell(nb, 1);
    matBlocks = cell(nb, 1);
    for ib = 1:nb
        datablockc = datablock{ib} ;
        if doubleCenter == 1
            bsize = size(datablock{ib}) ;
            d = bsize(1) ;
            n = bsize(2) ;
            datablockc = datablock{ib} - repmat(mean(datablock{ib},2), 1, n) - repmat(mean(datablock{ib},1), d, 1) + repmat(mean(mean(datablock{ib})), d, n) ; % BOTH CENTERING
        end
        [matBlocks{ib}, matLoadings{ib}] = MatReconstructMJ( ...
            datablockc, matJointV{ib}, matJointOrder{ib}, ...
            matJointRanks{ib});
    end
    
    outstruct.jointBasisMap = outMap;
    outstruct.matLoadings = matLoadings;
    outstruct.matBlocks = matBlocks;
    outstruct.keyIdxMap = keyIdxMap;
end


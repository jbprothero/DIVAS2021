function outstruct = DJIVEMainMJ_v1(datablock, paramstruct)
% DJIVEMainMJ_v1   DJIVE main functin, version 1
%
% Inputs:
%   datablock - cell array of d_k x n data matrices
%
%   paramstruct - a Matlab structure of input parameters (optional)
%                    Use: "help struct" and "help datatypes" to
%                         learn about these.
%                    Create one, using commands of the form:
%
%       paramstruct = struct('field1',values1,...
%                            'field2',values2,...
%                            'field3',values3) ;
%
%                          where any of the following can be used,
%                          these are optional, misspecified values
%                          revert to defaults
%
%    fields            values
%    dataname          nb x 1 cell array of string of data matrix's name
%    
%    nsim              number of bootstrap samples in Step 1
%
%   Copyright (c)  Meilei Jiang 2018

    % Initialize parameters
    nb = length(datablock); 
    dataname = cell(nb, 1);
    for ib = 1:nb
        dataname{ib} = ['datablock' num2str(ib)];
    end
    nsim = 400;
    theta0 = 45;
    optArgin = [];
    iprint = false;
    figdir = '';
    
    if exist('paramstruct', 'var')
        if isfield(paramstruct, 'dataname')
            dataname = paramstruct.dataname;
        end
        
        if isfield(paramstruct, 'nsim')
            nsim = paramstruct.nsim;
        end
        
        if isfield(paramstruct, 'theta0')
            theta0 = paramstruct.theta0;
        end
        
        if isfield(paramstruct, 'optArgin')
            optArgin = paramstruct.optArgin;
        end
        
        if isfield(paramstruct, 'iprint')
            iprint = paramstruct.iprint;
        end
        
        if isfield(paramstruct, 'figdir')
            figdir = paramstruct.figdir;
        end
    end
    
    
    % Step 1: Estimate signal space and perturbation angle
    [VBars, phiBars, ~, rBars] = ...
        DJIVESignalExtractJP(datablock, dataname, nsim, 0, 0, );

    % Step 2: Estimate joint ( and partially joint ) structure
    [outMap, keyIdxMap, ~, jointBlockOrder] = DJIVEJointStrucEstimateJP( ...
        VBars, phiBars, rBars, dataname, theta0, optArgin, iprint, figdir);

    % Step 3: Reconstruct DJIVE decomposition
    outstruct = DJIVEReconstructMJ(datablock, dataname, outMap, ...
        keyIdxMap, jointBlockOrder);

end
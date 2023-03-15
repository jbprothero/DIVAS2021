function outstruct = DJIVEMainJP(datablock, paramstruct, truth)
% DJIVEMainJP_v1   DJIVE main functin, version 1
% 
% Wrapper for three phases of DIVAS
% 
% DJIVESignalExtractJP chooses ranks and perturbation angles for each data
% block
%
% DJIVEJointStructureEstimateJP uses the ranks, subspaces, and angles from
% the previous step in an optimization problem to choose basis directions
% for each chunk of joint structure
%
% DJIVEReconstructMJ calculates corresponding loadings for each bundle of
% partially-shared joint structure
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
%    optArgin          Optimization paramters:
% {Tau_initial, Tau_max, Tau_multiplier, max_iterations,
% convergence_tolerance, slack_tolerance}
%
%    iprint            Boolean, whether optimization figures are produced
%
%    figdir            Folder location for optimization figures
%
%    colCent           Boolean, whether blocks should be object centered
%
%    rowCent           Boolean, whether blocks should be trait centered
% (as with other Marron software, rows are traits and columns are objects)
%
%    filterPerc        Percentage of RDB set as max for perturbation angle
%
%   Copyright (c)  Meilei Jiang 2018, Jack Prothero 2020

    disp('DIVAS Version 11-18-22 `Submission!`')
    % Initialize parameters
    nb = length(datablock); 
    dataname = cell(nb, 1);
    for ib = 1:nb
        dataname{ib} = ['datablock' num2str(ib)];
    end
    nsim = 400;
    theta0 = 45;
    optArgin = {0.5 1000 1.05 50 1e-3 1e-3};
    iprint = true;
    colCent = 0;
    rowCent = 0;
    figdir = '';
    filterPerc = 1 - (2/(1+sqrt(5))) ; % "Golden Ratio"
    %set a default
    noisepercentile = repmat(0.5, [1,nb]);
    
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
        
        if isfield(paramstruct, 'colCent')
            colCent = paramstruct.colCent;
        end
        
        if isfield(paramstruct, 'rowCent')
            rowCent = paramstruct.rowCent;
        end
        
        if isfield(paramstruct, 'filterPerc')
            filterPerc = paramstruct.filterPerc;
        end

        if isfield(paramstruct, 'noisepercentile')
            %vector with percentile of empiracle singular value to be
            %used in noise estimation specific to each block
            noisepercentile = paramstruct.noisepercentile;
        end

    end
    
    rowSpaces = cell(nb, 1);
    datablockc = cell(nb, 1);
    for ib = 1:nb
        rowSpaces{ib} = 0;
        datablockc{ib} = MatCenterJP(datablock{ib}, colCent, rowCent) ;
    end
    
    if exist('truth','var')
        for ib = 1:nb
            rowSpaces{ib} = orth(truth{ib}');
        end
    end
    
    % Step 1: Estimate signal space and perturbation angle
    [VBars, UBars, phiBars, psiBars, EHats, rBars, singVals, singValsHat, rSteps, VVHatCacheBars, UUHatCacheBars] = ...
        DJIVESignalExtractJP(datablockc, dataname, nsim, 0, colCent, rowCent, filterPerc, noisepercentile);
    
    delete(gcp('nocreate'))
    
    % Step 2: Estimate joint ( and partially joint ) structure
    [outMap, keyIdxMap, anglesMap, jointBlockOrder] = DJIVEJointStrucEstimateJPLoadInfo( ...
        VBars, UBars, phiBars, psiBars, rBars, datablockc, dataname, theta0, optArgin, iprint, figdir);
   
    % Step 3: Reconstruct DJIVE decomposition
    outstruct = DJIVEReconstructMJ(datablockc, dataname, outMap, ...
        keyIdxMap, jointBlockOrder, 0);
    
    outstruct.rBars = rBars ;
    outstruct.phiBars = phiBars ;
    outstruct.psiBars = psiBars ;
    outstruct.VBars = VBars ;
    outstruct.UBars = UBars ;
    outstruct.VVHatCacheBars = VVHatCacheBars ;
    outstruct.UUHatCacheBars = UUHatCacheBars ;
    outstruct.jointBasisMapRaw = outMap ;
    
    %{
    save SigExt.mat VBars UBars phiBars psiBars EHats rBars singVals singValsHat rSteps
    save cachedBootCull.mat VVHatCacheBars UUHatCacheBars
    save jointEstimate outMap keyIdxMap anglesMap jointBlockOrder
    save loadEstimate outstruct
    %}

end

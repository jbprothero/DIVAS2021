function [Vi, curRanks, angles] = BlockJointStrucEstimateJPLoadInfo(blockIn, datablock, dataname, ...
    VBars, UBars, phiBars, psiBars, rBars, curRanks, outMap, theta0, optArgin, iprint, ...
    figdir)
% BlockJointStrucEstimateMJ   Estimate a specific joint block basis
%   Detailed explanation goes here
%
% Inputs:
%   blockIn - a vector of logical, which indicating shared blocks
%   VBars - cell array of adjusted signal row spaces
%   phiBars - vector of perturbation angle for each data matrix
%   rBars - vector of adjusted signal ranks for each data matrix
%   curRanks - current cumulative ranks for each data block 
%   outMap - mapping between block index set and estimated partially shared
%            structure.
%   theta0 - angle between estimated spaces and optimized vector
%   optArgin - a cell array of optmization tuning parameters (optional): 
%              tau0, tau_max, mu, t_max, tol, delta
%
% Outputs:
%   Vi - estimated basis matrix
%   curRanks - updated cumulative ranks for each data blocks
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
    
    nb = length(blockIn);
    allIdx = 1:1:nb;
    blockIdx = allIdx(blockIn);
    blockName = dataname(blockIn);
    fprintf(strcat(strjoin({'Find joint structure shared only among', ...
        strjoin(blockName, ', ')}), '.\n'))
    
    n = size(VBars{1}, 1);
    blockLen = sum(blockIn);
    % set up optimization constrains
    Mo1 = [];
    Mo2 = [];
    Qc1 = cell(nb, 1);
    Qc2 = cell(nb, 1);
    Qc1Load = cell(nb, 1);
    Qc2Load = cell(nb, 1);
    precompute_load_matrix = cell(nb, 1);
    for ib = 1:nb
        d = size(UBars{ib}, 1) ;
        precompute_load_matrix{ib} = datablock{ib}' * UBars{ib} * UBars{ib}' * datablock{ib} ;
        if blockIn(ib)
            Mo2 = [Mo2, VBars{ib}];
            Qc1{ib} = eye(n, n);
            Qc2{ib} = VBars{ib} * VBars{ib}' / cosd(phiBars(ib))^2;
            Qc1Load{ib} = datablock{ib}' * datablock{ib};
            Qc2Load{ib} = precompute_load_matrix{ib} / cosd(psiBars(ib))^2;
        else
            % THIS PART SEEMS TO CONTROL THE PUSHING AWAY FROM NONINCLUDED
            % BLOCKS. IT WAS ORIGINALLY UNIMPLEMENTED BUT I ADDED IT
            % THINKING IT WAS MISTAKENLY LEFT OUT. FOR THE FUTURE, CONSIDER
            % ADDING AN ADJUSTABLE WEIGHT TO Mo1 TO CONTROL IMPORTANCE
            % Mo1 = [Mo1, VBars{ib}]; 
            Qc1{ib} = VBars{ib} * VBars{ib}' / cosd(phiBars(ib))^2;
            Qc2{ib} = eye(n, n);
            Qc1Load{ib} = zeros(n,n);
            Qc2Load{ib} = zeros(n,n);
            %{
            Qc1Load{ib} = precompute_load_matrix{ib} / cosd(psiBars(ib))^2;
            Qc2Load{ib} = datablock{ib}' * datablock{ib};
            %}
        end
    end
    if size(Mo1, 2) == 0
        Qo1 = 1e-6*ones(n, n);
    else
        Qo1 = Mo1 * Mo1';
    end
    if size(Mo2, 2) == 0
        Qo2 = 1e-6*ones(n, n);
    else
        Qo2 = Mo2 * Mo2';
    end
    
    
    % cross block constraints
    % CONSIDER EDITING THIS CHUNK TO FORCE 45 DEGREES BETWEEN ALL SUBSPACES
    % NOT JUST SUBSPACES IN DIFFERENT LAYERS
    % CURRENT BOUND BETWEEN ANY TWO SUBSPACES APPEARS TO BE MAX OF
    % PERTURBATION ANGLES AMONG BLOCKS INVOLVED IN EITHER SUBSPACE
    Vorth = [];
    Vnorth = [];
    for len = nb:-1:(blockLen+1)
        lenIdces = nchoosek(allIdx, len);
        nlen = size(lenIdces, 1);
        
        for i = 1:nlen
            bkIdx = lenIdces(i, :);
            bkIn = ismember(allIdx, bkIdx);
            t = Idx2numMJ(bkIn);
            if ~isKey(outMap, num2str(t))
                continue;
            end
            
            if all(ismember(blockIdx, bkIdx))
                Vorth = [Vorth, outMap(num2str(t))];
            else
                Vnorth = [Vnorth, outMap(num2str(t))];
            end
        end
    end

    if size(Vnorth, 1) > 0
        Qc1{end + 1} = Vnorth * Vnorth' / cosd(theta0)^2;
        Qc2{end + 1} = eye(n, n);
    end
    
    
    
    % starting optimization
    searchNext = true;
    
    if blockLen == 1
        proj = orth(Vorth);
        %left multiplication as projecting in column space
        [~,~,V0] = svds((datablock{blockIn}*(eye(size(Vorth,1))-proj*proj')), max(rBars(blockIn)));
    else
        [~,~,V0] = svds(Qo2 , max(rBars(blockIn)));
    end

    Vi = [];
    angles = 90*zeros(nb, size(V0,2));
    j = 0;
    while searchNext
        j = j + 1;
        if j > max(rBars(blockIn))
            break
            % THIS IS A HACKY FIX TO PREVENT A CRASH POSSIBILITY WHILE
            % ALLOWIING FOR DOUBLING UP BLOCK ENERGY FOR DIFFERENT BLOCK
            % COLLECTIONS. CONSIDER DOING SOMETHING SMARTER LONG-TERM
        end
        fprintf('Search Direction %d:\n', j)
        Vorth = [Vorth, Vi];
        output = penaltyCCPJPEarlyStopLoadInfo(V0(:,j), Qo1, Qo2, Qc1, Qc2, Qc1Load, Qc2Load, Vorth, optArgin);
        [opt_v, cache_v, ~, ~, converge] = output{:};
        if converge<=1
            angleHats = ccpOutAnalysisMJ(cache_v, VBars);
            %angleTrues = ccpOutAnalysisMJ(cache_v, rowSpaces);
            figname = strjoin({strjoin(blockName, '-'), '-joint-optV', ...
                num2str(j)}, '');
            ccpOutVisualMJ(angleHats, phiBars, dataname, iprint, figdir, figname);%, angleTrues);
            for ib = 1:nb
                angles(ib,j) = angleHats{ib}(end);
            end
        end
        if converge~=1
            fprintf('Direction %d does not converge. Stop searching current joint block.\n', j)
            break;
        end
        fprintf('Direction %d converges.\n', j)
        Vi = [Vi, opt_v];
        curRanks = curRanks + blockIn';
        %{
        if any(curRanks + blockIn' > rBars)
            fprintf('There is no room for searching next direction. Stop seraching current joint block.\n')
            searchNext = false;
        else
            fprintf('There is room for searching next direction. Continue...\n')
        end
        %}
    end
    
    angles = angles(:,1:(size(Vi,2)+1*~any(curRanks + blockIn' > rBars)));
    
end


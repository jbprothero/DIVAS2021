function [VBars, UBars, phiBars, psiBars, EHats, rBars, singVals, singValsHat, rSteps, VVHatCacheBars, UUHatCacheBars]= DJIVESignalExtractJP(datablock, ...
    dataname, nsim, iplot, colCent, rowCent, cull, noisepercentile, noiselvls)
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
    UBars = cell(nb, 1);
    EHats = cell(nb, 1);
    phiBars = 90 * ones(nb, 1);
    psiBars = 90 * ones(nb, 1);
    rBars = zeros(nb, 1);
    singVals = cell(nb, 1);
    singValsHat = cell(nb, 1);
    rSteps = cell(nb, 1);
    VVHatCacheBars = cell(nb, 1);
    UUHatCacheBars = cell(nb, 1);

    for ib = 1:nb
        fprintf('Signal estimation for %s\n', dataname{ib});
        datablockc = datablock{ib} ;
        d = size(datablockc, 1) ;
        n = size(datablockc, 2) ;
        percentile = noisepercentile(ib);
     
    
        if ~exist('noiselvls', 'var')
            [VBars{ib}, UBars{ib}, phiBars(ib), psiBars(ib), rBars(ib), EHats{ib}, singVals{ib}, singValsHat{ib}, rSteps{ib}, VVHatCacheBars{ib}, UUHatCacheBars{ib}] = ...
                MatSignalExtractJP(datablockc, dataname{ib}, nsim, colCent, rowCent, cull, percentile);
        else
            [VBars{ib}, UBars{ib}, phiBars(ib), psiBars(ib), rBars(ib), EHats{ib}, singVals{ib}, singValsHat{ib}, rSteps{ib}, VVHatCacheBars{ib}, UUHatCacheBars{ib}] = ...
                MatSignalExtractJP(datablockc, dataname{ib}, nsim, colCent, rowCent, cull, noiselvls{ib});
        end  
    end
    
    if iplot==1
        for ib = 1:nb
            singValsI = singVals{ib};
            singValsHatI = singValsHat{ib};
            rStepsI = rSteps{ib};
            rHat = sum(singValsHatI > 0);
            rBar = rBars(ib);
            matName = dataname{ib};
            figure
            topPlot = max(singValsI)+5;
            topRanges = topPlot*(1:length(rStepsI))/length(rStepsI);
            botRanges = topRanges - (topPlot/length(rStepsI));
            hold on
            plot(1:min(d,n), singValsI, '.b-', 'MarkerSize', 8)
            plot(1:min(d,n), singValsHatI, 'r-x')
            %plot([0 length(singVals)], [recovBound recovBound], 'g-')
            plot([rHat rHat], [0 topPlot], 'r-')
            plot([rBar rBar], [0 topPlot], 'k-')
            for i = 1:length(rStepsI)
                redProp = (length(rStepsI)-i)/(length(rStepsI)-1);
                plot([rStepsI(i) rStepsI(i)], [botRanges(i) topRanges(i)], '-', 'Color', [redProp 0 0])
            end
            xlim([0 length(singValsI)+1])
            ylim([0 topPlot])
            ylabel('Singular Value')
            title([matName ' Singular Value Shrinkage & Culling'])
        end
    end
end

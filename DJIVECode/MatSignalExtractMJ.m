function [VBar, phiBar, rBar, EHat] = MatSignalExtractMJ(X, matName, nsim)
% MatSignalExtractMJ   Matrix signal extraction
%   Estimate signal rank, signal row space, corresponding perturbation
%   angle and noise matrix. Adjust signals based on random direction angle.
%
% Inputs:
%   X - d x n data matrix
%   matName - string of matrix name
%   nsim - number of bootstrap samples
%
% Outputs:
%   VBar - adjusted signal row space
%   phiBar - adjusted perturbation angle
%   rBar - adjusted signal rank
%   EHat - estimated noise matrix
%
%   Copyright (c)  Meilei Jiang 2018

    singVals = svd(X);
    [d, n] = size(X);
    beta = min(n/d, d/n);
    [singValsHat,noiselvl] = optimal_shrinkage(singVals, beta, 'op');
    rHat = sum(singValsHat > 0);
    fprintf('Initial singal rank for %s is %d. \n', matName, rHat);
    recovBound = noiselvl*(sqrt(d)+sqrt(n));
    
    [UHat, ~, VHat] = svds(X, rHat);
    singValsTilde = singValsHat(1:rHat);
    AHat = UHat * diag(singValsTilde) * VHat';
    EHat = X - AHat; 
    
    randAngleCache = randDirAngleMJ(n, rHat, 1000);
    randAngle = quantile(randAngleCache, 0.05); % TRYING THIS AT 5 PERCENT
    
    % bootstrap estimation
    VVHatCache = cell(nsim, 1);
    PCAnglesCache = 90 * ones(nsim, rHat);
    
    parfor i = 1:nsim
        randV = orth(randn(n, rHat));
        randU = orth(randn(d, rHat));
        randX = randU * diag(singValsTilde) * randV' + EHat;
        [~, ~, randVHat] = svds(randX, rHat);
        VVHatCache{i} = randV' * randVHat;        
        PCAnglesCache(i, :) = acosd(min(abs(svd(VVHatCache{i})), 1));
    end
    
    % adjust signal rank by removing PCs with perturbation angle larger
    % than half value of random direction angle
    PCAngles = quantile(PCAnglesCache, 0.95, 1);
    validPC = (PCAngles < randAngle/2); % Ripe for changing
    rBar = sum(validPC > 0);
    fprintf('Adjusted singal rank for %s is %d. \n', matName, rBar);
    rSteps = [rHat rBar];
    % JP additions
    % Option 1: Repeated culling until every PC is below threshold
    %
    VVHatCacheBar = VVHatCache;
    while sum(validPC) < length(validPC)
        VVHatCacheBar = cell(nsim, 1);
        PCAnglesCacheBar = 90 * ones(nsim, rBar);
        singValsBar = singValsTilde(validPC);
        parfor i = 1:nsim
            randV = orth(randn(n, rBar));
            randU = orth(randn(d, rBar));
            randX = randU * diag(singValsBar) * randV' + EHat;
            [~, ~, randVHat] = svds(randX, rBar);
            VVHatCacheBar{i} = randV' * randVHat;        
            PCAnglesCacheBar(i, :) = acosd(min(svd(VVHatCacheBar{i}), 1));
        end
        PCAngles = quantile(PCAnglesCacheBar, 0.95, 1);
        validPC = (PCAngles < randAngle/2); % Ripe for changing
        rBar = sum(validPC > 0);
        rSteps = [rSteps rBar];
        fprintf('Adjusted singal rank for %s is %d. \n', matName, rBar);
    end
    maxAnglesBar = maxPerturbAngleMJ(ones(rBar,1), VVHatCacheBar);
    phiBar = quantile(maxAnglesBar, 0.95);
    rSteps(length(rSteps)) = [];
    %}
    
    % Option 2: Retain SVs associated with each remaining PC without
    % recalculating
    %{
    maxAnglesSequester = 90*ones(nsim,1) ;
    parfor i = 1:nsim
        tempSV = svd(VVHatCache{i}) ;
        maxAnglesSequester(i) = acosd(min(tempSV(validPC))) ;
    end
    phiBar = quantile(maxAnglesSequester, 0.95) ;
    
    %}
    % End of JP additions
    
    % Option 3: Old MJ Way, subset rows/columns of cached matrices and
    % recalculate svds of submatrices, doesn't appear to be doing the right
    % thing
    %{
    maxAngles = maxPerturbAngleMJ(validPC, VVHatCache);
    phiBar = quantile(maxAngles, 0.95);
    %}
    VBar = VHat(:, validPC);
    
    fprintf('Perturbation Angle for %s is %.1f. \n', matName, phiBar);
    
    % Singular Value Diagnostic Plot DON'T RUN THIS HERE
    %{
    figure
    topPlot = max(singVals)+5;
    topRanges = topPlot*(1:length(rSteps))/length(rSteps);
    botRanges = topRanges - (topPlot/length(rSteps));
    hold on
    plot(1:min(d,n), singVals, '.b-', 'MarkerSize', 8)
    plot(1:min(d,n), singValsHat, 'r-x')
    %plot([0 length(singVals)], [recovBound recovBound], 'g-')
    plot([rHat rHat], [0 topPlot], 'r-')
    plot([rBar rBar], [0 topPlot], 'k-')
    for i = 1:length(rSteps)
        redProp = (length(rSteps)-i)/(length(rSteps)-1);
        plot([rSteps(i) rSteps(i)], [botRanges(i) topRanges(i)], '-', 'Color', [redProp 0 0])
    end
    xlim([0 length(singVals)+1])
    ylim([0 topPlot])
    ylabel('Singular Value')
    title([matName ' Singular Value Shrinkage & Culling'])
    %}
end


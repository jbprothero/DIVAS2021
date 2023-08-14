function [VBar, UBar, phiBar, psiBar, rBar, EHat, singVals, singValsHat, rSteps, VVHatCacheBar, UUHatCacheBar] = MatSignalExtractJP(X, matName, nsim, colCent, rowCent, cull, percentile, noiselvl)
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

    
    
    [d, n] = size(X);
    mindn = min(d,n) ;
    % X = X - repmat(mean(X,2), 1, n) ; % Trying Column-Object Centering
    [UFull,singVals,VFull] = svd(X);
    singVals = diag(singVals) ;
    beta = min(n/d, d/n);
    if ~exist('cull', 'var')
        cull = 0.5 ;
    end
    
    
    if ~exist('noiselvl', 'var')
        [singValsHat,noiselvl] = optimal_shrinkage(singVals, beta, 'op', percentile);
    else
        if noiselvl == 'ks'
            noiselvl = ksOpt(singVals, beta) ;
        end
        singValsHat = optimal_shrinkage(singVals, beta, 'op', noiselvl);
    end
    rHat = sum(singValsHat > 0);
    fprintf('Initial singal rank for %s is %d. \n', matName, rHat);
    recovBound = noiselvl*(1+sqrt(beta));
    
    [UHat, ~, VHat] = svds(X, rHat);
    singValsTilde = singValsHat(1:rHat);
    AHat = UHat * diag(singValsTilde) * VHat';
    EHat = X - AHat ;
    EHatGood = UFull(:,(rHat+1):mindn) * diag(singVals((rHat+1):mindn)) * VFull(:,(rHat+1):mindn)' ;
    XRemaining = X - EHatGood ;

    % "Imputation" of missing energy
    imputedSingVals = zeros(1,rHat) ;
    for iter = 1:rHat
        perc = rand() ;
        marpas = PercentileMarcenkoPastur(beta, perc) ;
        imputedSingVals(iter) = sqrt(marpas);%*noiselvlShrinkage ;
    end
    EHatImpute = EHatGood + UHat * diag(imputedSingVals*noiselvl) * VHat' ;
    
    randAngleCache = randDirAngleMJ(n, rHat, 1000);
    randAngleCacheLoad = randDirAngleMJ(d, rHat, 1000);
    randAngle = quantile(randAngleCache, 0.05); % TRYING THIS AT 5 PERCENT
    randAngleLoad = quantile(randAngleCacheLoad, 0.05);
    
    % bootstrap estimation
    
    % adjust signal rank by removing PCs with perturbation angle larger
    % than half value of random direction angle
    
    rSteps = rHat;
    % JP additions
    % Option 1a: Repeated culling until every PC is below threshold
    %{
    PCAngles = quantile(PCAnglesCache, 0.95, 1);
    validPC = (PCAngles < randAngle/2); % Ripe for changing
    rBar = sum(validPC > 0);
    fprintf('Adjusted singal rank for %s is %d. \n', matName, rBar);
    rSteps = [rHat rBar];
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
    
    % Option 1b: Calculate maximum angle to bootstrap truth for every rank
    % between 1 and rHat of randVHat, choose rBar as the largest rank with
    % 95% below randAngle/2
    %
    PCAnglesCacheFullBoot = 90 * ones(nsim, rHat);
    PCAnglesCacheFullBootLoad = 90 * ones(nsim, rHat);
    
    fprintf('Progress Through Bootstraped Matrices:\n');
    fprintf(['\n' repmat('.',1,nsim) '\n\n']);
    
    parfor s = 1:nsim  
        randV = randn(n, rHat);
        if colCent
            randV = randV - repmat(mean(randV,1), n, 1);
        end
        randV = orth(randV);
        randU = randn(d, rHat);
        if rowCent
            randU = randU - repmat(mean(randU,1), d, 1);
        end
        randU = orth(randU);
        randX = randU * diag(singValsTilde) * randV' + EHatImpute;
        [randUHat, ~, randVHat] = svds(randX, rHat);
        for j=1:rHat
            PCAnglesCacheFullBoot(s,j) = acosd(min(svd(randV' * randVHat(:,1:j))));
            PCAnglesCacheFullBootLoad(s,j) = acosd(min(svd(randU' * randUHat(:,1:j)))) ;
        end
        fprintf('\b|\n');
    end
    
    rBar = sum(quantile(PCAnglesCacheFullBoot, 0.95, 1)<randAngle*cull) ;
    rBarLoad = sum(quantile(PCAnglesCacheFullBootLoad, 0.95, 1)<randAngleLoad*cull) ;
    %}
    
    % Option 1c: Calculate maximum angle to rank-varying bootstrap truth
    % for every rank between 1 and rHat, choose rBar as the largest rank
    % with 95% below randAngle/2
    %{
    PCAnglesCacheFullBoot = 90 * ones(nsim, rHat);
    
    fprintf('Progress Through Bootstraped Matrices:\n');
    
    parfor s = 1:nsim  
        randV = orth(randn(n, rHat));
        randU = orth(randn(d, rHat));
        randX = EHat ;
        for j=1:rHat
            randX = randX + singValsTilde(j) * randU(:,j) * randV(:,j)';
            [~, ~, randVHat] = svds(randX,j);
            PCAnglesCacheFullBoot(s,j) = acosd(min(svd(randV(:,1:j)' * randVHat(:,1:j))));
	    fprintf('%d\n',j);
        end
        fprintf('\b|s|\n');
    end
    %}
    fprintf(['Culled Rank is ' num2str(rBar) '.\n']) ;     
    validPC = quantile(PCAnglesCacheFullBoot, 0.95, 1)<randAngle*cull ;
    [~,minInd] = min(quantile(PCAnglesCacheFullBoot, 0.95, 1)) ;
    validPC(minInd) = 1 ;
    rBar = sum(validPC) ;
    phiBar = quantile(PCAnglesCacheFullBoot(:,rBar), 0.95) ;
    psiBar = quantile(PCAnglesCacheFullBootLoad(:,rBar), 0.95) ;
    
    VVHatCacheBar = cell(nsim, 1);
    UUHatCacheBar = cell(nsim, 1);
    singValsTildeBar = singValsTilde(1:rBar);
    parfor s = 1:nsim
        randV = randn(n, rBar);
        if colCent
            randV = randV - repmat(mean(randV,1), n, 1);
        end
        randV = orth(randV);
        randU = randn(d, rBar);
        if rowCent
            randU = randU - repmat(mean(randU,1), d, 1);
        end
        randU = orth(randU);
        randX = randU * diag(singValsTildeBar) * randV'+ EHat;
        [randUHat, ~, randVHat] = svds(randX, rBar);
        VVHatCacheBar{s} = randV' * randVHat;
        UUHatCacheBar{s} = randU' * randUHat;
    end
    
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
    UBar = UHat(:, validPC);
    
    fprintf('Perturbation Angle for %s is %.1f. \n', matName, phiBar);
    
    % Singular Value Diagnostic Plot
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


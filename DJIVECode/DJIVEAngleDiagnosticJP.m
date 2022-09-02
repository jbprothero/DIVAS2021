function DJIVEAngleDiagnosticJP(datablock, dataname, outstruct, randseed, titlestr)
%{
load('TCGA.mat')
load('tcgaSigExtBC04.mat')
load('jointEstimateBC04_noPushAway_CorrectedOpt.mat')
% load('jointReconstruct.mat')
load('cachedBootCullBC04.mat')
%}

rng(randseed)

outMap = outstruct.jointBasisMap ;
keyIdxMap = outstruct.keyIdxMap ;
rBars = outstruct.rBars ;
VBars = outstruct.VBars ;
UBars = outstruct.UBars ;
phiBars = outstruct.phiBars ;
psiBars = outstruct.psiBars ;
VVHatCacheBars = outstruct.VVHatCacheBars ;
UUHatCacheBars = outstruct.UUHatCacheBars ;


uniqueStr = titlestr ;
logENC = true ;
trueTraitKeys = keys(outstruct.jointBasisMap) ;

n = size(outMap(trueTraitKeys{1}),1);
nb = length(datablock);
ds = zeros(1,nb) ;
randAngleTraits = 90*ones(1,nb) ;
randAngleObjects = 90*ones(1,nb) ;
for ib = 1:nb
    trueObjKeys = keys(outstruct.matLoadings{ib}) ;
    [ds(ib),~] = size(outstruct.matLoadings{ib}(trueObjKeys{1})) ;
    randAngleCache = randDirAngleMJ(n, rBars(ib), 1000);
    randAngleCacheLoad = randDirAngleMJ(ds(ib), rBars(ib), 1000);
    randAngleTraits(ib) = quantile(randAngleCache, 0.05);
    randAngleObjects(ib) = quantile(randAngleCacheLoad, 0.05);
end


% CNS Plots
%{
cols = zeros(n,3) ;
markers = char(n,1) ;
for o = 1:n
    if strcmp('LumA',datasubtype(o))
        cols(o,:) = [0 0 1] ;
        markers(o) = '*';
    end
    if strcmp('LumB',datasubtype(o))
        cols(o,:) = [0 1 1] ;
        markers(o) = 'x';
    end
    if strcmp('Her2',datasubtype(o))
        cols(o,:) = [1 0 1] ;
        markers(o) = '+';
    end
    if strcmp('Basal',datasubtype(o))
        cols(o,:) = [1 0 0] ;
        markers(o) = '<';
    end
end
%}

%{
threeWay = outstruct.jointBasisMap('7') ;
paramstructScat = struct('npcadiradd',0, ...
                         'icolor', cols, ...
                         'isubpopkde', 1, ...
                         'markerstr', markers, ...
                         'datovlaymin', 0.3, ...
                         'datovlaymax', 0.8, ...
                         'legendcellstr', {{'LumA' 'LumB' 'Her2' 'Basl'}}, ...
                         'mlegendcolor', [0 0 1 ; 0 1 1 ; 1 0 1 ; 1 0 0], ...
                         'titlecellstr', {{'Joint Component Scatter Plots'}}, ...
                         'labelcellstr', {{'GE-CN-RPPA-Mut' 'GE-CN-RPPA-1' 'GE-CN-RPPA-2' 'GE-CN-RPPA-3' }}) ;
scatplotSM([outstruct.jointBasisMap('15') threeWay(:,1:3)]', eye(4), paramstructScat)
%}

% Subgroup Rank Diagram

% Identify each binary encoding by its number of blocks and its size as a
% subspace
rankIm = ones(nb,2^nb-1,3) ;
sspSize = zeros(1,2^nb-1) ;
sspCols = zeros(1,2^nb-1,3) ; 
numJ = zeros(1,2^nb-1) ;
sizeCols = [ 0 1 0 ; 1 0 0 ; 0.4 0.4 1 ; 0.5 0.5 0.5 ; 1 0 1 ] ;

for k = 1:(2^nb-1)
    key = num2str(k) ;
    if isKey(outMap, key)
        numJ(k) = length(keyIdxMap(key)) ;
        sspSize(k) = size(outMap(key),2) ;
        sspCols(:,k,:) = sizeCols(numJ(k),:) ;
        rankIm(keyIdxMap(key), k, :) =  repmat(sspCols(:,k,:),1,numJ(k),1);
    else
        bb = 1:nb ;
        binBlock = zeros(1,nb) ;
        rem = k ;
        for p = fliplr(1:nb)
            binBlock(p) = idivide(rem, int32(2^(p-1)), 'floor') ;
            rem = mod(rem, 2^(p-1)) ;
        end
        binBlock = logical(binBlock) ;
        kimk = bb(binBlock) ;
        keyIdxMap(key) = kimk ;
        numJ(k) = sum(binBlock) ;
        sspSize(k) = 0 ;
        sspCols(:,k,:) = [0.8 0.8 0.8] ;
        rankIm(kimk, k, :) =  repmat(sspCols(:,k,:),1,numJ(k),1);
    end
end

% Tweak the colors to weight them by subspace size
[numJ_s,szInd] = sort(numJ, 'descend') ;
sspSize_s = sspSize(szInd) ;
sspCols_s = sspCols(:,szInd,:) ;
rankIm_s_solid = rankIm(:,szInd,:) ;
rankIm_s = rankIm(:,szInd,:) ;
rankIm_s = [rankIm_s ones(nb,1,3)] ;

ncs = 1:(2^nb-1) ;
ins = 2^nb;
for b = fliplr(1:nb)
    ins = ins - nchoosek(nb,b) ;
    rankIm_s = [rankIm_s(:,1:ins,:) ones(nb,1,3) rankIm_s(:,(ins+1):size(rankIm_s,2),:)] ;
end
textInds = [] ;
runStart = 1 ;
runEnds = 1 ;
for b = 1:nb
    runEnds(b+1) = runEnds(b) + nchoosek(nb,b) ;
    textInds = [textInds runStart:(runStart-1+nchoosek(nb,b-1))] ;
    runStart = runStart+1+nchoosek(nb,b-1) ;
end
runEnds = runEnds(1:nb) ;

ncrs = zeros(1,nb);
midpts = zeros(1,nb);
xSum = 0;
xticklabs = cell(1,nb);
%rHats = zeros(1,nb);
for b = 1:nb
    ncrs(b) = nchoosek(nb,b-1);
    midpts(b) = mean(textInds((xSum+1):(xSum+ncrs(b)))) ;
    xSum = xSum+ncrs(b);
    xticklabs(nb-b+1) = cellstr([num2str(b) '-Way']);
    %rHats(b) = sum(singValsHat{b}>0) ;
end
xticklabs(nb+1) = cellstr('Totals');

% Plot!
%
figure
im = image(rankIm_s) ;
imax = ancestor(im, 'axes') ;
imyax = imax.YAxis ;
imyax.FontSize = 20;
title(['Rank Breakdown by Joint Structure: ' uniqueStr])
yticks(1:nb)
yticklabels(dataname)
xticks([midpts size(rankIm_s,2)])
xticklabels(xticklabs)
fontsize = 28 ;
hold on
for k = 1:(2^nb+nb)
    plot([k+0.5 k+0.5], [0 nb+1], 'k-')
end
xSum = 0 ;
for b = 1:nb
    for bb = 1:nb
        plot([textInds(xSum+1)-0.5 textInds(xSum+ncrs(b))+0.5], [bb-0.5 bb-0.5], 'k-')
    end
    text(2^nb+nb,b-0,{num2str(rBars(b))}, ...
        'FontSize',fontsize/2, 'HorizontalAlignment','center', 'Color','black')
    xSum = xSum+ncrs(b);
    plot([2^nb+nb-0.5 n^nb+nb+0.5], [b-0.5 b-0.5], 'k-')
end
for k = 1:(2^nb-1)
    k_s = szInd(k) ;
    for b = 1:nb
        if ismember(b, keyIdxMap(num2str(k_s)))
            text(textInds(k),b,num2str(sspSize(k_s)), 'HorizontalAlignment','center', 'FontSize',fontsize/2, 'Color','black');
        end
    end
end


%}

% CNS for every subset
%{
for i = 1:(2^nb-1)
    k = szInd(i) ;
    if ~isKey(outMap,num2str(k))
        continue
    end
    blocksInc = keyIdxMap(num2str(k)) ;
    blocksName = strjoin(dataname(blocksInc), '-') ;
    estSubspace = outMap(num2str(k)) ;
    if sspSize(k) > 1
        numPlots = min(4,sspSize(k)) ;
        firstTwo = estSubspace(:,1:numPlots) ;
        projmat = eye(numPlots) ;
        labelcellstr = {[blocksName '-V1'] [blocksName '-V2'] [blocksName '-V3'] [blocksName '-V4']};
        labelcellstr = labelcellstr(1:numPlots);
        paramstructScat = struct('npcadiradd',0, ...
                         'icolor', cols, ...
                         'isubpopkde', 1, ...
                         'markerstr', markers, ...
                         'datovlaymin', 0.3, ...
                         'datovlaymax', 0.8, ...
                         'legendcellstr', {{'LumA' 'LumB' 'Her2' 'Basl'}}, ...
                         'mlegendcolor', [0 0 1 ; 0 1 1 ; 1 0 1 ; 1 0 0], ...
                         'titlecellstr', {{[blocksName 'Joint Component Scatter Plots']}}, ...
                         'labelcellstr', {labelcellstr}) ;
        scatplotSM(firstTwo', projmat, paramstructScat)
    else
        first = estSubspace ;
        projmat = 1 ;
        labelcellstr = {{[blocksName '-V1']}} ;
        paramstructScat = struct('icolor', cols, ...
                         'isubpopkde', 1, ...
                         'markerstr', markers, ...
                         'datovlaymin', 0.3, ...
                         'datovlaymax', 0.8, ...
                         'legendcellstr', {{'LumA' 'LumB' 'Her2' 'Basl'}}, ...
                         'mlegendcolor', [0 0 1 ; 0 1 1 ; 1 0 1 ; 1 0 0], ...
                         'titlestr', {{[blocksName ' Joint Component Scatter Plots']}}, ...
                         'xlabelstr', labelcellstr) ;
        projplot1SM(first', projmat, paramstructScat)
    end
    
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [ 0 0 14 8 ] ;
    print([blocksName '-BC-noPushAway-cns'], '-dpng', '-r0')
    close
    
end
%}

% Finding bounds for chosen directions

%
%
thetaTwoStars = containers.Map ;
upperBounds = containers.Map ;
lowerBounds = containers.Map ;
midPoints = containers.Map ;

thetaTwoStarsLoad = containers.Map ;
upperBoundsLoad = containers.Map ;
lowerBoundsLoad = containers.Map ;
midPointsLoad = containers.Map ;

for k = 1:(2^nb-1)
    fprintf('Starting binary code %d\n', k);
    if isKey(outMap, num2str(k))
        jointVec = outMap(num2str(k)) ;
        rs = size(jointVec,2) ;
        thetaTwos = zeros(nb, rs) ;
        thetaTwosLoad = zeros(nb, rs) ;
        midpoints = zeros(nb, rs) ;
        midpointsLoad = zeros(nb, rs) ;
        for ib = 1:nb
            if sum(sum(VVHatCacheBars{ib}{1} ~= 0)) ~= 0
                nsim = size(VVHatCacheBars{ib}, 1) ;
                %effSingVals = singValsHat{ib}(1:rBars(ib));
                
                if isKey(outstruct.matLoadings{ib}, num2str(k))
                    loadVec = outstruct.matLoadings{ib}(num2str(k)) ;
                    loadVec = takeNormOfEachColumnJP(loadVec) ;
                else
                    loadVec = zeros(size(datablock{ib},1),rs) ;
                end
                
                omegaHat = VBars{ib}' * jointVec;
                omegaHatLoad = UBars{ib}' * loadVec;
                
                % min/max to keep in domain
                for r = 1:rs
                    midpoints(ib, r) = acosd(min(1,max(-1,svds(omegaHat(:,r), 1)))) ;
                    midpointsLoad(ib, r) = acosd(min(1,max(-1,svds(omegaHatLoad(:,r), 1)))) ;
                end
                
                thetaTwoStarsBoot = zeros(nsim,rs);
                precos = zeros(nsim,rs);
                thetaTwoStarsBootLoad = zeros(nsim,rs);
                precosLoad = zeros(nsim,rs);
                
                bootCache = VVHatCacheBars{ib};
                bootCacheLoad = UUHatCacheBars{ib};
                
                fprintf('Progress Through Bootstraped Matrices:\n');
                fprintf(['\n' repmat('.',1,nsim) '\n\n']);
                parfor i = 1:nsim
                    bootMat = bootCache{i};
                    bootMatLoad = bootCacheLoad{i};
                    
                    precos(i,:) = sqrt(sum((bootMat*omegaHat).^2,1))./sqrt(sum(omegaHat.^2,1));
                    precosLoad(i,:) = sqrt(sum((bootMatLoad*omegaHatLoad).^2,1))./sqrt(sum(omegaHatLoad.^2,1));
                    
                    thetaTwoStarsBoot(i,:) = acosd(sqrt(sum((bootMat*omegaHat).^2,1))./sqrt(sum(omegaHat.^2,1)));
                    thetaTwoStarsBootLoad(i,:) = acosd(sqrt(sum((bootMatLoad*omegaHatLoad).^2,1))./sqrt(sum(omegaHatLoad.^2,1)));
                    
                    fprintf('\b|\n');
                end
                thetaTwos(ib,:) = quantile(thetaTwoStarsBoot, 0.95, 1);
                thetaTwosLoad(ib,:) = quantile(thetaTwoStarsBootLoad, 0.95, 1);
            end
        end
        thetaTwoStars(num2str(k)) = thetaTwos ;
        %midpoints = anglesMap(num2str(k)) ;
        upperbounds = thetaTwos + midpoints(:,1:rs) ;
        upperbounds(upperbounds > 90) = 90 ;
        upperBounds(num2str(k)) = upperbounds ;
        lowerbounds = midpoints(:,1:rs) - repmat(phiBars, 1, rs) ;
        lowerbounds(lowerbounds < 0) = 0 ;
        lowerBounds(num2str(k)) = lowerbounds ;
        midPoints(num2str(k)) = midpoints ;
        
        thetaTwoStarsLoad(num2str(k)) = thetaTwosLoad ;
        upperboundsLoad = thetaTwosLoad + midpointsLoad ;
        upperboundsLoad(upperboundsLoad > 90) = 90 ;
        upperBoundsLoad(num2str(k)) = upperboundsLoad ;
        lowerboundsLoad = midpointsLoad - repmat(psiBars, 1, rs) ;
        lowerboundsLoad(lowerboundsLoad < 0) = 0 ;
        lowerBoundsLoad(num2str(k)) = lowerboundsLoad ;
        midPointsLoad(num2str(k)) = midpointsLoad ;
    end
    fprintf('Finished binary code %d\n', k) ;
end

%
% Plotting lower & upper bounds for scores

lastRow = ones(1, size(rankIm_s,2), 3) ;
lastRow(:,textInds,:) = sspCols_s ;

figure
im = image([rankIm_s ; lastRow]) ;
imax = ancestor(im, 'axes') ;
imyax = imax.YAxis ;
imyax.FontSize = 20;
title(['Joint Structure Score Diagnostics: ' uniqueStr])
yticks(1:nb)
yticklabels(dataname)
text(0.48,nb+1,{'Effective';'Number';'of Cases'}, 'HorizontalAlignment','right')
xticks([midpts 2^nb+nb])
xticklabels(xticklabs)
fontsize = 28 ;
sepBar = 0.5 ;
hold on
for k = 1:(2^nb+nb)
    plot([k+0.5 k+0.5], [0 nb+2], 'k-')
end
xSum = 0 ;
for b = 1:nb
    for bb = 1:nb
        plot([textInds(xSum+1)-0.5 textInds(xSum+ncrs(b))+0.5], [bb-0.5 bb-0.5], 'k-')
    end
    text(2^nb+nb,b-0.4,{num2str(rBars(b))}, ...
        'FontSize',fontsize/2, 'HorizontalAlignment','center', 'Color','black')
    xSum = xSum+ncrs(b);
    plot([2^nb+nb-0.5 n^nb+nb+0.5], [b-0.5 b-0.5], 'k-')
end
plot([-0.5 size(lastRow,2)+0.5], [nb+0.5 nb+0.5], '-', 'Color',[0.4 0.4 0.4], 'LineWidth',1)
plot([-0.5 size(lastRow,2)+0.5], [nb+0.6 nb+0.6], '-', 'Color',[0.4 0.4 0.4], 'LineWidth',1)
plot([-0.5 size(lastRow,2)+0.5], [nb+0.55 nb+0.55], '-', 'Color',[0.4 0.4 0.4], 'LineWidth',11)
for k = 1:(2^nb-1)
    k_s = szInd(k) ;
    k_p = textInds(k) ;
    if isKey(outMap, num2str(k_s))
        ddExtend = 0.5 ;
        if ismember(k,runEnds)
            for b = 1:nb
                text(k_p+0.55, -phiBars(b)/90+b+0.5, [num2str(round(phiBars(b),1)) char(176)], 'FontSize',fontsize/4)
                text(k_p+0.75, -randAngleTraits(b)./90+b+0.5, [num2str(round(randAngleTraits(b),1)) char(176)], 'FontSize',fontsize/4)
            end
            ddExtend = 0.72 ;
        end
        xspan = linspace(k_p-0.5, k_p+0.5, sspSize(k_s)+2) ; 
        xspan = xspan(2:(sspSize(k_s)+1)) ;
        angleMat = midPoints(num2str(k_s)) ;
        if size(angleMat,2) > sspSize(k_s)
            angleMat(:,size(angleMat,2)) = [] ;
        end
        upperMat = upperBounds(num2str(k_s)) ;
        lowerMat = lowerBounds(num2str(k_s)) ;
        thetaTwos = thetaTwoStars(num2str(k_s)) ;
        for b = 1:nb
            plot(xspan, -angleMat(b,:)/90+b+0.5, 'kx', 'MarkerSize',4)
            if sum(sum(VVHatCacheBars{b}{1} ~= 0)) ~= 0
                plot(xspan, -upperMat(b,:)/90+b+0.5, 'k.', 'MarkerSize',9)
                %plot(xspan, -lowerMat(b,:)/90+b+0.5, 'w.-')
            end
            plot([k_p-0.5 k_p+0.5], -repelem(phiBars(b),2)/90+b+0.5, 'k--')
            plot([k_p-0.5 k_p+ddExtend], -repelem(randAngleTraits(b),2)/90+b+0.5, 'k-.')
            %plot([k_s-0.5 k_s+0.5], [b-0.5 b-0.5], 'k-')
        end
    end
    for b = 1:nb
        if ismember(b, keyIdxMap(num2str(k_s)))
            text(k_p,b-0.4,num2str(sspSize(k_s)), 'HorizontalAlignment','center', 'FontSize',fontsize/2, 'Color','black');
        end
    end
end

%
essRoot = ceil(n^(1/4)) ;
for k=1:(2^nb-1)
    k_s = szInd(k) ;
    k_p = textInds(k) ;
    if isKey(outMap, num2str(k_s))
        xspan = linspace(k_p-0.5, k_p+0.5, sspSize(k_s)+2) ; 
        xspan = xspan(2:(sspSize(k_s)+1)) ;
        effSampSizes = 1./sum(outMap(num2str(k_s)).^4) ;
        
        if logENC == true
            plot(xspan, -.9*log(effSampSizes)/log(n)+nb+1+0.5, 'k+', 'MarkerSize',4)
            plot([k_p-0.5 k_p+0.5], -.9*repelem(log(essRoot),2)/log(n)+nb+1+0.5, 'k--')
            plot([k_p-0.5 k_p+0.5], -.9*repelem(log(essRoot^2),2)/log(n)+nb+1+0.5, 'k--')
            plot([k_p-0.5 k_p+0.5], -.9*repelem(log(essRoot^3),2)/log(n)+nb+1+0.5, 'k--')
            % plot([k_p-0.5 k_p+0.5], -repelem(log(256),2)/log(n)+b+0.5, 'k--')
            if ismember(k,runEnds)
                text(k_p+0.55, -.9*log(essRoot)/log(n)+nb+1+0.5, num2str(essRoot), 'FontSize',fontsize/4)
                text(k_p+0.55, -.9*log(essRoot^2)/log(n)+nb+1+0.5, num2str(essRoot^2), 'FontSize',fontsize/4)
                text(k_p+0.55, -.9*log(essRoot^3)/log(n)+nb+1+0.5, num2str(essRoot^3), 'FontSize',fontsize/4)
            end
        else
            plot(xspan, -.9*effSampSizes/n+nb+1+0.5, 'k+', 'MarkerSize',4)
            plot([k_p-0.5 k_p+0.5], repelem(-.9*0.75+nb+1+0.5,2), 'k--')
            plot([k_p-0.5 k_p+0.5], repelem(-.9*0.50+nb+1+0.5,2), 'k--')
            plot([k_p-0.5 k_p+0.5], repelem(-.9*0.25+nb+1+0.5,2), 'k--')
            % plot([k_p-0.5 k_p+0.5], -repelem(log(256),2)/log(n)+b+0.5, 'k--')
            if ismember(k,runEnds)
                text(k_p+0.55, -.9*0.75+nb+1+0.5, num2str(round(0.75*n)), 'FontSize',fontsize/4)
                text(k_p+0.55, -.9*0.50+nb+1+0.5, num2str(round(0.50*n)), 'FontSize',fontsize/4)
                text(k_p+0.55, -.9*0.25+nb+1+0.5, num2str(round(0.25*n)), 'FontSize',fontsize/4)
            end
        end
    end
    %{
    for b = 1:nb
        if ismember(b, keyIdxMap(num2str(k_s)))
            text(k_p,nb+1+0.4,num2str(sspSize(k_s)), 'HorizontalAlignment','center', 'FontSize',fontsize/2, 'Color','black');
        end
    end
    %}
end

%}

%
% Single Cell for Demonstration
%{
sampleCol = 0.4*ones(1,1,3) ;
sampleCol(:,:,3) = 1 ;
figure
im = image(sampleCol) ;
imax = ancestor(im, 'axes') ;
imyax = imax.YAxis ;
imyax.FontSize = 20;
yticks(1)
yticklabels('X1')
xticks(1)
xticklabels('3-Way')
fontsize = 28 ;
angleMat = anglesMap('7') ;
upperMat = upperBounds('7') ;

hold on
plot(1, -angleMat(1,1)/90+1+0.5, 'kx', 'MarkerSize',8)
plot(1, -upperMat(1,1)/90+1+0.5, 'k.', 'MarkerSize',18)
plot([0.5 1.5], -repelem(phiBars(1),2)/90+1+0.5, 'k--', 'LineWidth',2)
plot([0.5 1+ddExtend], -repelem(randAngleTraits(1),2)/90+1+0.5, 'k-.', 'LineWidth',2)
text(1+0.51, -phiBars(1)/90+1+0.5, [num2str(round(phiBars(1),1)) char(176)], 'FontSize',fontsize/2)
text(1+0.51, -randAngleTraits(1)./90+1+0.5, [num2str(round(randAngleTraits(1),1)) char(176)], 'FontSize',fontsize/2)
%text(1,1-0.4,'1', 'HorizontalAlignment','center', 'FontSize',fontsize/2, 'Color','black');
%}


% Plotting lower & upper bounds for loadings

figure
im = image([rankIm_s ; lastRow]) ;
imax = ancestor(im, 'axes') ;
imyax = imax.YAxis ;
imyax.FontSize = 20;
title(['Joint Structure Loadings Diagnostics: ' uniqueStr])
yticks(1:nb)
yticklabels(dataname)
text(0.48,nb+1,{'Effective';'Contribution';'of Traits'}, 'HorizontalAlignment','right')
xticks([midpts 2^nb+nb])
xticklabels(xticklabs)
fontsize = 28 ;
sepBar = 0.5 ;
hold on
for k = 1:(2^nb+nb)
    plot([k+0.5 k+0.5], [0 nb+2], 'k-')
end
xSum = 0 ;
for b = 1:nb
    for bb = 1:nb
        plot([textInds(xSum+1)-0.5 textInds(xSum+ncrs(b))+0.5], [bb-0.5 bb-0.5], 'k-')
    end
    text(2^nb+nb,b-0.4,{num2str(rBars(b))}, ...
        'FontSize',fontsize/2, 'HorizontalAlignment','center', 'Color','black')
    xSum = xSum+ncrs(b);
    plot([2^nb+nb-0.5 n^nb+nb+0.5], [b-0.5 b-0.5], 'k-')
end
plot([-0.5 size(lastRow,2)+0.5], [nb+0.5 nb+0.5], '-', 'Color',[0.4 0.4 0.4], 'LineWidth',1)
plot([-0.5 size(lastRow,2)+0.5], [nb+0.6 nb+0.6], '-', 'Color',[0.4 0.4 0.4], 'LineWidth',1)
plot([-0.5 size(lastRow,2)+0.5], [nb+0.55 nb+0.55], '-', 'Color',[0.4 0.4 0.4], 'LineWidth',11)
for k = 1:(2^nb-1)
    k_s = szInd(k) ;
    k_p = textInds(k) ;
    if isKey(outMap, num2str(k_s))
        ddExtend = 0.5 ;
        if ismember(k,runEnds)
            for b = 1:nb
                text(k_p+0.55, -psiBars(b)/90+b+0.5, [num2str(round(psiBars(b),1)) char(176)], 'FontSize',fontsize/4)
                text(k_p+0.75, -randAngleObjects(b)./90+b+0.5, [num2str(round(randAngleObjects(b),1)) char(176)], 'FontSize',fontsize/4)
            end
            ddExtend = 0.72 ;
        end
        xspan = linspace(k_p-0.5, k_p+0.5, sspSize(k_s)+2) ; 
        xspan = xspan(2:(sspSize(k_s)+1)) ;
        angleMat = midPointsLoad(num2str(k_s)) ;
        if size(angleMat,2) > sspSize(k_s)
            angleMat(:,size(angleMat,2)) = [] ;
        end
        upperMat = upperBoundsLoad(num2str(k_s)) ;
        lowerMat = lowerBoundsLoad(num2str(k_s)) ;
        thetaTwos = thetaTwoStarsLoad(num2str(k_s)) ;
        for b = 1:nb
            if sum(sum(UUHatCacheBars{b}{1} ~= 0)) ~= 0
                plot(xspan, -upperMat(b,:)/90+b+0.5, 'k.', 'MarkerSize',9)
                %plot(xspan, -lowerMat(b,:)/90+b+0.5, 'w.-')
            end
            if isKey(outstruct.matLoadings{b}, num2str(k_s))
                plot(xspan, -angleMat(b,:)/90+b+0.5, 'kx', 'MarkerSize',4)
                plot([k_p-0.5 k_p+0.5], -repelem(psiBars(b),2)/90+b+0.5, 'k--')
                plot([k_p-0.5 k_p+ddExtend], -repelem(randAngleObjects(b),2)/90+b+0.5, 'k-.')
                %plot([k_s-0.5 k_s+0.5], [b-0.5 b-0.5], 'k-')
            end
        end
    end
    for b = 1:nb
        if ismember(b, keyIdxMap(num2str(k_s)))
            text(k_p,b-0.4,num2str(sspSize(k_s)), 'HorizontalAlignment','center', 'FontSize',fontsize/2, 'Color','black');
        end
    end
end

% Loadings Effective Traits
for k=1:(2^nb-1)
    k_s = szInd(k) ;
    k_p = textInds(k) ;
    if isKey(outMap, num2str(k_s))
        xspan = linspace(k_p-0.5, k_p+0.5, sspSize(k_s)+2) ; 
        xspan = xspan(2:(sspSize(k_s)+1)) ;
        for ib = 1:nb
            if isKey(outstruct.matLoadings{ib}, num2str(k_s))
                loadVec = outstruct.matLoadings{ib}(num2str(k_s)) ;
                loadVec = takeNormOfEachColumnJP(loadVec) ;
                effSampSizes = 1./sum(loadVec.^4) ;
                text(xspan, -.9*effSampSizes/ds(ib)+nb+1+0.5, num2str(ib), 'FontSize',fontsize/4)
            end
        end
        
        plot([k_p-0.5 k_p+0.5], repelem(-.9*0.75+nb+1+0.5,2), 'k--')
        plot([k_p-0.5 k_p+0.5], repelem(-.9*0.50+nb+1+0.5,2), 'k--')
        plot([k_p-0.5 k_p+0.5], repelem(-.9*0.25+nb+1+0.5,2), 'k--')
        % plot([k_p-0.5 k_p+0.5], -repelem(log(256),2)/log(n)+b+0.5, 'k--')
        if ismember(k,runEnds)
            text(k_p+0.55, -.9*0.75+nb+1+0.5, '75%', 'FontSize',fontsize/4)
            text(k_p+0.55, -.9*0.50+nb+1+0.5, '50%', 'FontSize',fontsize/4)
            text(k_p+0.55, -.9*0.25+nb+1+0.5, '25%', 'FontSize',fontsize/4)
        end
    end
    %{
    for b = 1:nb
        if ismember(b, keyIdxMap(num2str(k_s)))
            text(k_p,nb+1+0.4,num2str(sspSize(k_s)), 'HorizontalAlignment','center', 'FontSize',fontsize/2, 'Color','black');
        end
    end
    %}
end

% Loadings Energy
%{
for k=1:(2^nb-1)
    k_s = szInd(k) ;
    k_p = textInds(k) ;
    if isKey(outMap, num2str(k_s))
        xspan = linspace(k_p-0.5, k_p+0.5, sspSize(k_s)+2) ; 
        xspan = xspan(2:(sspSize(k_s)+1)) ;
        loadEnergies = zeros(nb,sspSize(k_s)) ;
        for ib = 1:nb
            if isKey(outstruct.matLoadings{ib}, num2str(k_s))
                loadVec = outstruct.matLoadings{ib}(num2str(k_s)) ;
                loadVec = takeNormOfEachColumnJP(loadVec) ;
                for u = 1:sspSize(k_s)
                    loadEnergies(ib,u) = sum(sum((loadVec(:,u) * loadVec(:,u)' * datablockc{ib}).^2)) ;
                end
                loadEnergies(ib,:) = loadEnergies(ib,:)/sum(sum(datablockc{ib}.^2)) ;
                text(xspan, -.9*loadEnergies(ib,:)+nb+1+0.5, num2str(ib), 'FontSize',fontsize/4)
            end
        end
        
        plot([k_p-0.5 k_p+0.5], repelem(-.9*0.75+nb+1+0.5,2), 'k--')
        plot([k_p-0.5 k_p+0.5], repelem(-.9*0.50+nb+1+0.5,2), 'k--')
        plot([k_p-0.5 k_p+0.5], repelem(-.9*0.25+nb+1+0.5,2), 'k--')
        % plot([k_p-0.5 k_p+0.5], -repelem(log(256),2)/log(n)+b+0.5, 'k--')
        if ismember(k,runEnds)
            text(k_p+0.55, -.9*0.75+nb+1+0.5, '75%', 'FontSize',fontsize/4)
            text(k_p+0.55, -.9*0.50+nb+1+0.5, '50%', 'FontSize',fontsize/4)
            text(k_p+0.55, -.9*0.25+nb+1+0.5, '25%', 'FontSize',fontsize/4)
        end
    end
    %{
    for b = 1:nb
        if ismember(b, keyIdxMap(num2str(k_s)))
            text(k_p,nb+1+0.4,num2str(sspSize(k_s)), 'HorizontalAlignment','center', 'FontSize',fontsize/2, 'Color','black');
        end
    end
    %}
end
%}

% Plotting "Effective Sample Size" for each direction
%{
figure
im = image([rankIm_s]) ;
imax = ancestor(im, 'axes') ;
imyax = imax.YAxis ;
imyax.FontSize = 20;
title(['TCGA Joint Structure Effective Sample Sizes: ' uniqueStr])
yticks(1:nb)
yticklabels(dataname)
xticks([midpts 2^nb+nb])
xticklabels(xticklabs)
fontsize = 28 ;
sepBar = 0.5 ;
hold on
for k = 1:(2^nb+nb)
    plot([k+0.5 k+0.5], [0 nb+1], 'k-')
end
xSum = 0 ;
for b = 1:nb
    for bb = 1:nb
        plot([textInds(xSum+1)-0.5 textInds(xSum+ncrs(b))+0.5], [bb-0.5 bb-0.5], 'k-')
    end
    text(2^nb+nb,b+0.4,{num2str(rBars(b))}, ...
        'FontSize',fontsize/2, 'HorizontalAlignment','center', 'Color','black')
    xSum = xSum+ncrs(b);
    plot([2^nb+nb-0.5 n^nb+nb+0.5], [b-0.5 b-0.5], 'k-')
end
for k=1:(2^nb-1)
    k_s = szInd(k) ;
    k_p = textInds(k) ;
    if isKey(outMap, num2str(k_s))
        xspan = linspace(k_p-0.5, k_p+0.5, sspSize(k_s)+2) ; 
        xspan = xspan(2:(sspSize(k_s)+1)) ;
        effSampSizes = 1./sum(outMap(num2str(k_s)).^4) ;
        for b = 1:nb
            plot(xspan, -log(effSampSizes)/log(n)+b+0.5, 'k+', 'MarkerSize',4)
            plot([k_p-0.5 k_p+0.5], -repelem(log(5),2)/log(n)+b+0.5, 'k--')
            plot([k_p-0.5 k_p+0.5], -repelem(log(25),2)/log(n)+b+0.5, 'k--')
            plot([k_p-0.5 k_p+0.5], -repelem(log(125),2)/log(n)+b+0.5, 'k--')
            % plot([k_p-0.5 k_p+0.5], -repelem(log(256),2)/log(n)+b+0.5, 'k--')
            if b==1 && k==1
                text(k_p+0.55, -log(5)/log(n)+b+0.5, num2str(5), 'FontSize',fontsize/4)
                text(k_p+0.55, -log(25)/log(n)+b+0.5, num2str(25), 'FontSize',fontsize/4)
                text(k_p+0.55, -log(125)/log(n)+b+0.5, num2str(125), 'FontSize',fontsize/4)
            end
        end
    end
    for b = 1:nb
        if ismember(b, keyIdxMap(num2str(k_s)))
            text(k_p,b+0.4,num2str(sspSize(k_s)), 'HorizontalAlignment','center', 'FontSize',fontsize/2, 'Color','black');
        end
    end
end
%}

end

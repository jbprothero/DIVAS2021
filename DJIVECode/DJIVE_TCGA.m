function DJIVE_TCGA(cullingPercentage)

cull = 0.5 ;
if nargin == 1
    cull = cullingPercentage ;
end

cullstr = num2str(cull) ;
cullstr(cullstr=='.') = [] ;

load('TCGA.mat') ;
parpool([12 20]) ;

% Step 1: Estimate Signal Subspaces and Cache Bootstrap Matrices
[VBars, UBars, phiBars, EHats, rBars, singVals, singValsHat, rSteps, VVHatCacheBar, UUHatCacheBar] = ...
        DJIVESignalExtractJP(datablock, dataname, 400, 0, 1, 1, cull);
    
save(['tcgaSigExtBC' cullstr '.mat'], 'VBars', 'UBars', 'phiBars', 'EHats', 'rBars', 'singVals', 'singValsHat', 'rSteps')
save(['cachedBootCullBC' cullstr '.mat'], 'VVHatCacheBar', 'UUHatCacheBar')

% Step 2: Estimate joint ( and partially joint ) structure
[outMap, keyIdxMap, anglesMap, jointBlockOrder] = DJIVEJointStrucEstimateJP( ...
    VBars, phiBars, rBars, dataname, 45, {0.5 1000 1.03 100 1e-3 1e-3}, 1, '\Figures');

save(['jointEstimateBC' cullstr '_noPushAway_CorrectedOpt.mat'], 'outMap', 'keyIdxMap', 'anglesMap', 'jointBlockOrder')

% Step 3: Reconstruct DJIVE decomposition
outstruct = DJIVEReconstructMJ(datablock, dataname, outMap, ...
    keyIdxMap, jointBlockOrder);

save(['jointReconstructBC' cullstr '_CorrectedOpt.mat'], 'outstruct')

end

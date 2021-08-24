function sigma = ksOpt(singVals, betaShrinkage)

    sigmaMin = median(singVals) / (1 + sqrt(betaShrinkage)) ;
    sigmaMax = 2 * max(singVals) / (1 + sqrt(betaShrinkage)) ;
    numGrid = 200 ;
    
    cands = sigmaMin:((sigmaMax-sigmaMin)/(numGrid-1)):sigmaMax ;
    objVals = zeros(1,numGrid) ;
    
    for ii = 1:numGrid
        sigmaCand = cands(ii) ;
        noiseSingVals = singVals(singVals < sigmaCand * (1 + betaShrinkage)) ;
        card = length(noiseSingVals) ;
        
        absVals = zeros(1,card) ;
        
        for jj = 1:card
            absVals(jj) = abs(incMarPas((noiseSingVals(jj)/sigmaCand)^2, betaShrinkage, 0) - (jj - 0.5)/card) ;
        end
        
        objVals(ii) = max(absVals) + 1/(2*card) ;
        disp(['Finished noise estimation candidate ' num2str(ii)])
    end

    [~,minInd] = min(objVals) ;
    sigma = cands(minInd) ;
end

function I = incMarPas(x0,beta,gamma)

%{
    if beta > 1
        error('betaBeyond');
    end
%}
    topSpec = (1 + sqrt(beta))^2;
    botSpec = (1 - sqrt(beta))^2;
    MarPas = @(x) IfElse((topSpec-x).*(x-botSpec) >0, ...
                         sqrt((topSpec-x).*(x-botSpec))./(beta.* x)./(2 .* pi), ...
                         0);
    if gamma ~= 0
       fun = @(x) (x.^gamma .* MarPas(x));
    else
       fun = @(x) MarPas(x);
    end
    I = integral(fun,x0,topSpec);
    
    function y=IfElse(Q,point,counterPoint)
        y = point;
        if any(~Q)
            if length(counterPoint) == 1
                counterPoint = ones(size(Q)).*counterPoint;
            end
            y(~Q) = counterPoint(~Q);
        end
        
    end
end
function [singValsHat, noiselvl, rHat] = SingValShrinkPlotJP(X, paramstruct)

% SingValShrinkPlotJP
% Input a data matrix X
% Output a vector of its shrunken singular values according to 
% optimal_shrinkage from Gavish & Donoho.
% matName is an optional argument to input the name of X for plotting &
% console printing.
% numShow is the number of singular values to include in the plot. 
% (default 2*rHat)
% iPlot=1 will create a plot of the empirical & shrunken singular values
% (default 1)
% iLog=1 will plot & return the natural logarithm of singular values 
% instead 
% iLog=2 will plot squared singular values
% (default 0)
% iSave=1 will automatically save a copy of the plotted figure to the
% current directory as a .png file. (default 0)
% noiseLvl: User supplied noise level. If not provided the default is the
% median singular value divided by the square root of the marcenko pastur
% median
% numShow: Number of singular values to show on the output plot. Default is
% 2 times the shrunken rank
% iThresh: Type of thresholding: 0 no thresholding
%                                1 hard thresholding
%                                2 soft thresholding
%                                3 (default) G-D optimal operator norm loss
%                                thresholding

matName = 'Data Matrix' ;
iPlot = 1 ;
iLog = 0 ;
iSave = 0 ;
noiselvlFlag = 0 ;
numShowFlag = 0 ;
iThresh = 3 ;

if nargin > 1
    if isfield(paramstruct,'matName')     
        matName = getfield(paramstruct,'matName') ; 
    end 
    if isfield(paramstruct,'iPlot')     
        iPlot = getfield(paramstruct,'iPlot') ; 
    end
    if isfield(paramstruct,'iSave')     
        iSave = getfield(paramstruct,'iSave') ; 
    end
    if isfield(paramstruct,'noiselvl')     
        noiselvl = getfield(paramstruct,'noiselvl') ; 
        noiselvlFlag = 1 ;
    end
    if isfield(paramstruct,'numShow')     
        numShow = getfield(paramstruct,'numShow') ;
        numShowFlag = 1 ;
    end
    if isfield(paramstruct,'iThresh')     
        iThresh = getfield(paramstruct,'iThresh') ; 
    end
    if isfield(paramstruct,'iLog')     
        iLog = getfield(paramstruct,'iLog') ; 
    end
end

singVals = svd(X);
[d, n] = size(X);
beta = min(n/d, d/n);
hardThresh = sqrt(2*(beta+1)+(8*beta)/(beta+1+sqrt(beta^2+14*beta+1))) ;
softThresh = 1+sqrt(beta) ;

if ~noiselvlFlag
    [singValsHat,noiselvl] = optimal_shrinkage(singVals, beta, 'op');
    if iThresh == 1
        singValsHat = (singVals).*((singVals/noiselvl)>hardThresh) ;
    end
    if iThresh == 2
        singValsHat = noiselvl*max(((singVals/noiselvl)-softThresh), 0) ;
    end 
else
    if iThresh == 1
        singValsHat = (singVals).*((singVals/noiselvl)>hardThresh) ;
    end
    if iThresh == 2
        singValsHat = noiselvl*max(((singVals/noiselvl)-softThresh), 0) ;
    end 
    if iThresh == 3
        singValsHat = optimal_shrinkage(singVals, beta, 'op', noiselvl);
    end
end
rHat = sum(singValsHat > 0);
fprintf('Initial signal rank for %s is %d. \n', matName, rHat);
recovBound = noiselvl*beta^(1/4);
empirBound = noiselvl*(1+ sqrt(beta));
if iThresh == 1
    empirBound = noiselvl*hardThresh ;
end

[UHat, ~, VHat] = svds(X, rHat);
singValsTilde = singValsHat(1:rHat);
AHat = UHat * diag(singValsTilde) * VHat';
EHat = X - AHat;

singValsE = svd(EHat) ;
numQuants = length(singValsE) ;
quants = zeros(numQuants,1) ;
for q = 1:numQuants
    quants(q) = PercentileMarcenkoPastur(beta,q/(numQuants+1)) ;
end

%{
nSim = 400 ;
invSimulNums = rand(nSim,numQuants) ;
sampleMarPas = zeros(nSim,numQuants) ;
parfor s = 1:nSim
    for q = 1:numQuants
        sampleMarPas(s,q) = PercentileMarcenkoPastur(beta,invSimulNums(s,q)) ;
    end
    fprintf('Finished Row %d \n', s)
end
sampleMarPas = sort(sampleMarPas,2) ;
sampleMarPas = sampleMarPas - repmat(mean(sampleMarPas,2),1,numQuants) ;
sampleMarPas = sampleMarPas ./ repmat(std(sampleMarPas,[],2),1,numQuants) ;
lowerEnvelope = quantile(sampleMarPas, 0.025, 1) ;
upperEnvelope = quantile(sampleMarPas, 0.975, 1) ;
%}

%{
nSim = 100 ;
sampleMarPas = zeros(nSim,numQuants) ;
parfor s = 1:nSim
    tempSingVals = svd(normrnd(0,1,d,n)).^2 ;
    sampleMarPas(s,:) = (sort(tempSingVals)-mean(tempSingVals))/std(tempSingVals) ;
    fprintf('Finished Row %d \n', s)
end
%}

if ~numShowFlag
    numShow = min(2*rHat, length(singVals));
end

ylab = "Singular Value";
if iThresh == 1
    main = 'Singular Value Hard Thresholding';
end
if iThresh == 2
    main = 'Singular Value Soft Thresholding';
end
if iThresh == 3
    main = 'Singular Value Shrinkage';
end
if iLog==1
    singVals = log10(singVals);
    singValsHat = log10(singValsHat);
    recovBound = log10(recovBound);
    empirBound = log10(empirBound);
    ylab = "Log Singular Value";
    main = ['Log10 ' main];
end
if iLog==2
    singVals = singVals.^2;
    singValsHat = singValsHat.^2;
    recovBound = recovBound^2;
    empirBound = empirBound^2;
    ylab = "Eigenvalue";
    main = 'Eigenvalue Shrinkage';
end

if iPlot==1
    
    hold on
    topPlot = max(singVals)+1;
    botPlot = min(max([min(singVals); -5]), -0.5);
    plot(1:min(d,n), singVals, '.b-', 'MarkerSize', 8) % Empirical singvals
    plot(1:min(d,n), singValsHat, 'r-x') % Shrunken singvals
    plot([0 min(d,n)+1], [empirBound empirBound], 'r--') % empirical nonzero min
    plot([0 min(d,n)+1], [recovBound recovBound], 'g-') % theoretical recoverable min
    xlim([0 numShow+1])
    ylim([botPlot topPlot])
    ylabel(ylab)
    title([matName ' ' main])
    
    if iSave==1
        fig = gcf;
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [ 0 0 14 8 ] ;
        print([matName 'singValShrink'], '-dpng', '-r0')
        close
    end
    
    %{
    figure 
    hold on
    stdQuants = (sort(quants)-mean(quants))/std(quants) ;
    stdSingValsE = (sort(singValsE)-mean(singValsE))/std(singValsE) ;
    sqLim = ceil(max([abs(stdQuants) ; abs(stdSingValsE)])) ;
    plot(stdQuants, stdSingValsE, '.m-', 'MarkerSize', 8) % Error singvals
    plot(stdQuants, lowerEnvelope, 'b-')
    plot(stdQuants, upperEnvelope, 'b-')
    plot([-sqLim sqLim], [-sqLim sqLim], 'k--')
    plot([-sqLim sqLim], [0 0], 'k-')
    plot([0 0], [-sqLim sqLim], 'k-')
    xlabel('Theoretical Quantiles')
    ylabel('Error Singular Values')
    title([matName ' Noise Diagnostic'])
    %}
end

end

function out = PercentileMarcenkoPastur(beta, perc)
    MarPas = @(x) 1-incMarPas(x,beta,0);
    lobnd = (1 - sqrt(beta))^2;
    hibnd = (1 + sqrt(beta))^2;
    change = 1;
    while change && (hibnd - lobnd > .001)
      change = 0;
      x = linspace(lobnd,hibnd,5);
      for i=1:length(x)
          y(i) = MarPas(x(i));
      end
      if any(y < perc)
         lobnd = max(x(y < perc));
         change = 1;
      end
      if any(y > perc)
         hibnd = min(x(y > perc));
         change = 1;
      end
    end
    out = (hibnd+lobnd)./2;
end

function I = incMarPas(x0,beta,gamma)
    if beta > 1
        error('betaBeyond');
    end
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


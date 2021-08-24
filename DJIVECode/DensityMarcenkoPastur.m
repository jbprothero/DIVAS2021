function densityPoints = DensityMarcenkoPastur(beta, iSave)

betastr = num2str(beta) ;
betasavestr = betastr ;
betasavestr(regexp(betastr, '[.]')) = [] ;

topSpec = (1 + sqrt(beta))^2;
botSpec = (1 - sqrt(beta))^2;
points = linspace(botSpec, topSpec, 1000) ;
MarPas = @(x) IfElse((topSpec-x).*(x-botSpec) >0, ...
    sqrt((topSpec-x).*(x-botSpec))./(beta.* x)./(2 .* pi), ...
    0);

densityPoints = MarPas(points) ;
figure
plot([botSpec points topSpec], [0 densityPoints 0], '-r')
xlabel("Sample Covariance Matrix Eigenvalue")
ylabel("Density")
title(['Marchenko-Pastur Density for \beta = ' betastr])

if beta > 1
    figure
    subplot(1,2,1)
    plot([botSpec points topSpec], [0 densityPoints 0], '-r')
    xlabel("Sample Covariance Matrix Eigenvalue")
    ylabel("Density")
    title(['Marchenko-Pastur Density for \beta = ' betastr])
    hold on
    plot(0, 1-(1/beta), '.r', 'MarkerSize',10)
    plot([0 0], [0 1-(1/beta)], '-r')
    subplot(1,2,2)
    plot([botSpec points topSpec], [0 densityPoints 0], '-r')
    xlabel("Sample Covariance Matrix Eigenvalue")
    ylabel("Density")
    title(['Non-Zero Component Marchenko-Pastur Density for \beta = ' betastr])
else
    figure
    plot([botSpec points topSpec], [0 densityPoints 0], '-r')
    xlabel("Sample Covariance Matrix Eigenvalue")
    ylabel("Density")
    title(['Marchenko-Pastur Density for \beta = ' betastr])
end

if iSave == 1
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [ 0 0 14 8 ] ;
    print(['MarPasDensity' betasavestr], '-dpng', '-r0')
    close
end

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
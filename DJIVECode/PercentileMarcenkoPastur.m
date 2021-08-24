function out = PercentileMarcenkoPastur(beta, perc)

    nonzeroArea = 1 ;
    if beta >= 1
        nonzeroArea = 1/beta ;
        perc = perc*nonzeroArea ;
    end
    MarPas = @(x) nonzeroArea-incMarPas(x,beta,0);
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
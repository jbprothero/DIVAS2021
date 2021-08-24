function output = penaltyCCPJPEarlyStop(v0, Qo1, Qo2, Qc1, Qc2, Vorth, optArgin)
% penaltyCCP   Penalty convex concave procedure algorithm
%   Detailed explanation goes here
%
% Inputs:
%   v0 - Initial direction
%   Qo1 - a matrix of the 1st quadratic form in the objective function
%   Qo2 - a matrix of the 2nd quadratic form in the objective function
%   Qc1 - a nc x 1 cell array of matrices of the 1st quadratic form in the 
%         contraints
%   Qc2 - a nc x 1 cell array of matrices of the 2nd quadratic form in the 
%         contraints
%   Vorth - a matrix of orthogonal constraints
%   optArgin - a cell array of optmization tunning parameters: 
%              tau0, tau_max, mu, t_max, tol, delta
%
% Outputs:
%   outputArg1 - Description
%   outputArg2 - Description
%
%   Copyright (c)  Meilei Jiang 2018

    numvarargs = length(optArgin);
    optargs = {0.5 1000 1.05 200 1e-3 1e-3};
    optargs(1:numvarargs) = optArgin;
    [tau0, tau_max, mu, t_max, tol, delta] = optargs{:};
    
    if size(Vorth, 1) == 0
        Vorth = zeros(size(Qo1, 1), 1);
    end
    
    cache_v = cell(t_max, 1);
    cache_cvx_objval = inf*ones(t_max, 1);
    cache_slack = cell(t_max, 1);
    
    cache_v{1} = v0;
    nc = length(Qc1);
    cache_slack{1} = inf*ones(nc + 2, 1);
    
    converge = 0; % flag to ensure two steps without change for convergence
    tau = tau0;
    
    updatePrintFormat = ['%d ' repmat('%f ', 1, length(Qc1)+3) '\n'] ;
    
    for t = 2:t_max
        if mod(t, 10) == 0
            disp(['Iteration ' num2str(t)])
        end
        result = ccpSubOptJPEarlyStop(cache_v{t-1}, Qo1, Qo2, Qc1, Qc2, Vorth, tau);
        [cache_v{t}, cache_slack{t}, cache_cvx_objval(t)] = result{:};
        if any(isnan(cache_v{t}))
            disp(['NaN solution appears in iteration' num2str(t)])
            t = t-1;
            break
        end
        cache_v{t} = cache_v{t}/norm(cache_v{t});
        curr_objval = cache_cvx_objval(t) + tau * sum(cache_slack{t});
        pre_objval = cache_cvx_objval(t-1) + (tau/mu) * sum(cache_slack{t-1}); % TYPO????? _objval(t) instead of _objval(t-1)
        fprintf(updatePrintFormat, t, cache_slack{t}, curr_objval-pre_objval)
        curr_slack = cache_slack{t} ;
        if ((norm(curr_objval-pre_objval) < tol) && ...
                (sum(cache_slack{t}) <= delta))
            if converge == 1
                break
            else
                converge = 1;
            end
        else
            converge = 0;
        end
        tau = min(mu*tau, tau_max);
    end
    cache_v = cache_v(1:t);
    cache_cvx_objval = cache_cvx_objval(1:t);
    cache_slack = cache_slack(1:t);
    opt_v = cache_v{end};
    output = {opt_v, cache_v, cache_cvx_objval, cache_slack, converge};
end


function result = ccpSubOptJPEarlyStop(v0, Qo1, Qo2, Qc1, Qc2, Vo, tau)
% ccpSubOptMJ   Sub-optmization problem in penalty CCP algorithm
%   quadratic programming.
%
% Inputs:
%   v0 - Initial direction
%   Qo1 - a matrix of the 1st quadratic form in the objective function
%   Qo2 - a matrix of the 2nd quadratic form in the objective function
%   Qc1 - a nc x 1 cell array of matrices of the 1st quadratic form in the 
%         contraints
%   Qc2 - a nc x 1 cell array of matrices of the 2nd quadratic form in the 
%         contraints
%   Vo  - a matrix of orthogonal constraints
%   tau - multiplier of slack variables
%
% Outputs:
%   result - tuple of optimized direction, slack variable value, objective
%   value
%
%   Copyright (c)  Meilei Jiang 2018

    nc = length(Qc1);
    n = length(v0);
    [~, ro] = size(Vo);
    cvx_begin quiet;
        variables v(n) slack(nc+2)
        minimize(quad_form(v, Qo1) - 2*v0'*Qo2*v + quad_form(v0, Qo2) ...
            + tau*sum(slack))
        subject to
            for ic = 1:nc
                quad_form(v, Qc1{ic}) - 2*v0'*Qc2{ic} * v + quad_form(v0, Qc2{ic}) <= slack(ic);
            end
             v'*v - 1 <= slack(nc + 1);
            1 - 2*v0'*v + v'*v <= slack(nc + 2);
            slack >= 0;
            Vo' * v == zeros(ro, 1);            
    cvx_end
    cvx_objval = v'*Qo1*v - 2*v0'*Qo2*v + v0'*Qo2*v0;
    result = {v, slack, cvx_objval};
end


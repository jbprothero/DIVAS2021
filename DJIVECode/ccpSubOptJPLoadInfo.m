function result = ccpSubOptJPLoadInfo(v0, Qo1, Qo2, Qc1, Qc2, Qc1Load, Qc2Load, Vo, tau)
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
    ncload = length(Qc1Load);
    n = length(v0);
    d = ones(1,ncload);
    v0load = cell(ncload,1);
    loadslackscale = ones(1,ncload);
    for ipre = 1:ncload
        if sum(sum(Qc1Load{ipre} ~= 0))
            loadslackscale(ipre) = svds(Qc1Load{ipre},1) ;
        end
    end
    
    %{
    for i = 1:ncload
        d(i) = size(datablock{i},1) ;
        v0load{i} = datablock{i} * v0 ;
    end
    dlim = max(d) ;
    %}
    [~, ro] = size(Vo);
    cvx_begin quiet;
        variables v(n) slack(nc + ncload + 2)
        minimize(quad_form(v, Qo1) - 2*v0'*Qo2*v + quad_form(v0, Qo2) ...
            + tau*sum(slack))
        subject to
            for ic = 1:nc
                quad_form(v, Qc1{ic}) - 2*v0'*Qc2{ic} * v + quad_form(v0, Qc2{ic}) <= slack(ic); %no scaling here because Qc1/Qc2 have SV<=1
            end
            for ic = 1:ncload
                quad_form(v, Qc1Load{ic}) - 2*v0'*Qc2Load{ic}*v + quad_form(v0, Qc2Load{ic}) <= slack(nc+ic)/loadslackscale(ic); %scaling here because QcLoad1 (X'UU'X) has big SVs
            end
            v'*v - 1 <= slack(nc + ncload + 1);
            1 - 2*v0'*v + v'*v <= slack(nc + ncload + 2); %trying upweighting norm constraints to prevent counterintuitive solution paths.
            slack >= 0;
            Vo' * v == zeros(ro, 1);            
    cvx_end
    cvx_objval = v'*Qo1*v - 2*v0'*Qo2*v + v0'*Qo2*v0;
    result = {v, slack, cvx_optval, cvx_status};
end


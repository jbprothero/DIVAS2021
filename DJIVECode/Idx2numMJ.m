function t = Idx2numMJ(blockIn)
% Idx2numMJ   Convert logical index array to number (binary to decimal)
%   Take blockIN as binary representation and convert it to demcimal
%   number.
%
% Inputs:
%   blockIn - a vector of logical to indicate which block idx is in
%
% Outputs:
%   t - converted number
%
%   Copyright (c)  Meilei Jiang 2018

    nb = length(blockIn);
    t=0;
    for i = 1:nb
        if blockIn(i)
            t = t + 2^(i-1);
        end
    end    
end


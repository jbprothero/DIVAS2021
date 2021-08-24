function angle = projAngleMJ(x, V)
% projAngleMJ   Calculate the angle between vector x and V
%   Detailed explanation goes here
%
% Inputs:
%   x - a n x 1 vector
%   V - a n x r basis matrix
%
% Outputs:
%   angle - projected angle of x on V
%
%   Copyright (c)  Meilei Jiang 2018

    px = V * V' * x;
    angle = acosd(px' * x / (norm(px) * norm(x)) );
end


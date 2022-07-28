%
% Copyright (c) 2015, Mostapha Kalami Heris & Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "LICENSE" file for license terms.
%
% Project Code: YPEA101
% Project Title: Implementation of Binary Genetic Algorithm in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Cite as:
% Mostapha Kalami Heris, Binary and Real-Coded Genetic Algorithms in MATLAB (URL: https://yarpiz.com/23/ypea101-genetic-algorithms), Yarpiz, 2015.
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function [y1, y2] = SinglePointCrossover(x1, x2)

    nVar = numel(x1);
    
    c = randi([1 nVar-1]);
    
    y1 = [x1(1:c) x2(c+1:end)];
    y2 = [x2(1:c) x1(c+1:end)];

end
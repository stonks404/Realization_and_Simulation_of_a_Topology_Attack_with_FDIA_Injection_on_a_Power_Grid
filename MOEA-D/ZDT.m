%
% Copyright (c) 2015, Mostapha Kalami Heris & Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "LICENSE" file for license terms.
%
% Project Code: YPEA124
% Project Title: Implementation of MOEA/D
% Muti-Objective Evolutionary Algorithm based on Decomposition
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Cite as:
% Mostapha Kalami Heris, MOEA/D in MATLAB (URL: https://yarpiz.com/95/ypea124-moead), Yarpiz, 2015.
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function z = ZDT(x)

    n = numel(x);

    f1 = x(1);
    
    g = 1+9/(n-1)*sum(x(2:end));
    
    h = 1-sqrt(f1/g);
    
    f2 = g*h;
    
    z = [f1
       f2];

end
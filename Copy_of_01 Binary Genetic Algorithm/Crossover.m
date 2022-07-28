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

function [y1, y2] = Crossover(x1, x2)

    pSinglePoint = 0.1;
    pDoublePoint = 0.2;
    pUniform = 1-pSinglePoint-pDoublePoint;
    
    METHOD = RouletteWheelSelection([pSinglePoint pDoublePoint pUniform]);
    
    switch METHOD
        case 1
            [y1, y2] = SinglePointCrossover(x1, x2);
            
        case 2
            [y1, y2] = DoublePointCrossover(x1, x2);
            
        case 3
            [y1, y2] = UniformCrossover(x1, x2);
            
    end

end
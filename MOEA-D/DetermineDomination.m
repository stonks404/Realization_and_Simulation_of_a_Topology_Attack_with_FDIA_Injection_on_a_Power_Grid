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

function pop = DetermineDomination(pop)

    nPop = numel(pop);

    for i = 1:nPop
        pop(i).IsDominated = false;
    end
    
    for i = 1:nPop
        for j = i+1:nPop
            if Dominates(pop(i), pop(j))
                pop(j).IsDominated = true;
                
            elseif Dominates(pop(j), pop(i))
                pop(i).IsDominated = true;
                
            end
        end
    end

end
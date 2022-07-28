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

function y = Crossover(x1, x2, params)

    alpha = randi([0 1], size(x1));
    
    y = alpha.*x1+(1-alpha).*x2;
    y2 = alpha.*x2+(1-alpha).*x1;
    
    
    display('FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF');

%     gamma = params.gamma;
%     VarMin = params.VarMin;
%     VarMax = params.VarMax;
%     
%     alpha = unifrnd(-gamma, 1+gamma, size(x1));
%     
%     y = alpha.*x1+(1-alpha).*x2;
% 
%     y = round(y,0);
% 
%     %y = min(max(y, VarMin), VarMax);
% 
%     cond = 1;
%     %v = randi([0 1], VarSize);
%     while cond
%         gamma = params.gamma;
%         VarMin = params.VarMin;
%         VarMax = params.VarMax;
%     
%         alpha = unifrnd(-gamma, 1+gamma, size(x1));
%     
%         y = alpha.*x1+(1-alpha).*x2;
% 
%         y = round(y,0);
%         if (sum(y(1:33)) < 4)
%             cond = 0;
%         end
%     end
% display(x1);
% display(x2);
%       y = x1 & x2;
% 
%       num_modified_lines_omega_r = length(find(y(1:33,1) == 1));
%       num_modified_lines_omega_a = length(find(y(34:35,1) == 1)); 
%         
% 
%       index = find(y == 0);
%       index_len = length(index);
%       n = randi(index_len);
% 
%       y(n) = 1;
% 
%     display(y);


    
end
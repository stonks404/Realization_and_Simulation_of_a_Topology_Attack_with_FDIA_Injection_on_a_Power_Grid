function [newSwarm, adjustSuccess] = adjustSwarm(popSize, swarm, funcName_adjustInd, userObj)

adjustSuccess = ones(1, popSize); %flags to indicate whether the countries can satisfy the constraints;

for i=1:popSize
    cockroach = swarm(i, :);
    
    str_adjustInd = [funcName_adjustInd, '(cockroach, userObj)'];
    result = eval(str_adjustInd);
    if(result == -1) %if adjust constraint not sccess;
        swarm(i, :) = cockroach;
        adjustSuccess(1, i) = 0;
    else
        swarm(i, :) = result;
        adjustSuccess(1, i) = 1;
    end
end

newSwarm = swarm;

end
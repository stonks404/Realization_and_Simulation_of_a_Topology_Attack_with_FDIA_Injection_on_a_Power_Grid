function [newCockroach, fitness] = evaluateCockroach(cockroach, bounds, types, funcName_adjustInd, funcName_fitness, userObj)

maxFitness = 999999999999999999999999999999;
indToEval = cockroach;

dim = length(indToEval);

for i=1:dim
   value = cockroach(1, i);
   type = types(1, i);
   if(type==1)  %if the datatype of this dimension is binary;
       middle = (bounds(1, i) + bounds(2, i))/2;
       if(value<middle)
           value = 0;
       else
           value = 1;
       end
   elseif (type==2) %if the datatype of this dimension is integer;
       newValue = round(value);
       if(newValue < bounds(1, i))
            newValue = ceil(bounds(1, i));
       elseif (newValue > bounds(2, i))
            newValue = floor(bounds(2, i));
       end
       value = newValue;
   end
   
   indToEval(1, i) = value;
end

str_adjustInd = [funcName_adjustInd, '(indToEval, userObj)'];
[newInd, newUserObj] = eval(str_adjustInd);
if(newInd == -1)
    fitness = maxFitness;
else
    userObj = newUserObj;
    indToEval = newInd;
    str_fitness = [funcName_fitness, '(indToEval, userObj)'];
    fitness = eval(str_fitness);
end

newCockroach = indToEval;



end
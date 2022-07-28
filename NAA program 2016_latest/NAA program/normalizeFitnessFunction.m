%calculate the normalized fitness values of the shelter leaders;
%Inputs: sortedLeaderFitnesses-(1*shelterNum) vector. The sorted fitness values of the leaders;
%        fitness-the fitness value of a certain shelter leader;
%Output: value-the normalized fitness value of a certain shelter leader;
%Author: Fengji Luo
%Date: 07/2015
function value = normalizeFitnessFunction(sortedLeaderFitnesses, fitness, base)

n = length(sortedLeaderFitnesses);
distance = fitness - base;


sum_distance = 0;
for i=1:n
    tmp_distance = sortedLeaderFitnesses(i) - base;
    sum_distance = sum_distance + tmp_distance;
end

if(sum_distance==0)
    value = 1/n;
else
    value = 1- (distance/sum_distance);
end

end
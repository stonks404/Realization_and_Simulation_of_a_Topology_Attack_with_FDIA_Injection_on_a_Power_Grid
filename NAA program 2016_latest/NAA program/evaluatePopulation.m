%evaluate the fitness values of the whole population.
%Inputs: popSize-Integer. The number of the individuals;
%        swarm-(popSize*dimension) matrix. Rrepresenting the individuals;
%        bounds-(2*dimension) matrix. The boundaries of each
%                dimension;
%        types-(1*dimension) vector. The data types of each dimension:
%        funcName_adjustInd-The name of the individual adjustment function .m file;
%        funcName_fitness-String. The name of the fitness function .m file;
%        userObj-Structure. Contain some extra information (if needed)
%                which are needed when calculate the fitness value; 
%Output: fitnesses-(1*popSize) vector. Entry (1, i) represents the
%        calcualted fitness value of the ith individual;
%Author: Fengji Luo
%Date:   07/2015
function [newSwarm, fitnesses] = evaluatePopulation(popSize, swarm, bounds, types, funcName_adjustInd, funcName_fitness, userObj)

fitnesses = zeros(1, popSize);  %the fitness values of the population;
newSwarm = swarm;
for i=1:popSize  %firstly, evaluate the population;
    ind = swarm(i, :);
    [newCockroach, fitness] = evaluateCockroach(ind, bounds, types, funcName_adjustInd, funcName_fitness, userObj);
    newSwarm(i, :) = newCockroach;
    fitnesses(1, i) = fitness;
end
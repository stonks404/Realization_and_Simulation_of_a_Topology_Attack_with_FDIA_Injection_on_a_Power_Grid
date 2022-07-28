%initialize the shelters;
%Inputs: swarm-(popSize*dimension) matrix. The whole population;
%        dimension-Integer. The dimension of each individual;
%        bounds-(2*dimension) matrix. The boundaries of each
%                dimension;
%        types-(1*dimension) vector. The data types of each dimension:
%        popSize-Integer. Number of the individuals;
%        shelterNum-Integer. Number of the shelters;
%        shelterCaps-(1*shelterNum) vector. The capacity of the shelters;
%        funcName_adjustInd-The name of the individual adjustment function .m file;
%        funcName_fitness-The name of the fitness function .m file;
%        userObj-Structure. Contain some extra information (if needed)
%                which are needed when calculate the fitness value; 
%Output: swarmFitnesses-(1*popSize). Fitness values of the individuals;
%        shelterIndexes-(1*popSize). The shelter indexes of the
%                        individuals;
%        shelterLeaders-(shelterNum, dimension) matrix. Each row represents
%                        a shelter leader vector;
%        normalizedLeaderFitnesses-(1*shelterNum). Normalized fitness
%                                   values of the shelter leaders;
%Author: Fengji Luo
%Date: 07/2015
function [newSwarm, swarmFitnesses, shelterIndexes, shelterLeaders, normalizedLeaderFitnesses] = formShelters(swarm, dimension, bounds, types, popSize, shelterNum, shelterCaps, funcName_adjustInd, funcName_fitness, userObj)

shelterIndexes = zeros(1, popSize);  %0-corner, >0-the index of the belonged sub population;
for i=1:popSize
    shelterIndexes(1,i) = -1; %initialize the shelter index of each individual as -1;
end


shelterLeaders = zeros(shelterNum, dimension);
[newSwarm, swarmFitnesses] = evaluatePopulation(popSize, swarm, bounds, types, funcName_adjustInd, funcName_fitness, userObj); %evaluate the whole population;
[values, indexes] = sort(swarmFitnesses);  %sort the fitness values of the individuals;

%set the shelter leaders;
for i=1:shelterNum
    originIndex = indexes(i);
    shelterLeaders(i, :) = newSwarm(originIndex, :);
    shelterIndexes(1,originIndex) = i;
end

%set the shelter leader fitness values;
shelterFitness = zeros(1, shelterNum);
for i=1:shelterNum
    shelterFitness(1, i) = values(i);
end

%set the normalized shelter leader fitness values;
normalizedLeaderFitnesses = zeros(1, shelterNum);
for i=1:shelterNum
    normalizedLeaderFitnesses(1,i) = normalizeFitnessFunction(shelterFitness, shelterFitness(i), values(shelterNum+1));
end


%assign the individuals to the shelters;
start = 1;
for i=1:shelterNum
    cap = shelterCaps(i);  %capacity of the shelter i;
    count = 0;
    for j=start:popSize
        if(shelterIndexes(j)==-1)
            shelterIndexes(j) = i;
            count = count+1;
        end
        if(count >= cap)
            start = j;
            break;
        end
    end
end

end
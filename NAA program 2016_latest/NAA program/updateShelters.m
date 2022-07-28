%Update the shelters;;
%Inputs: dimension-Integer. Dimension of each individual;
%        swarm-(popSize*dimension) matrix. Set of the individuals;
%        shelterNum-Integer. Number of the shelters;
%        shelterIndexes-(1*popSize) vector. Shelter indexes of the
%                        individuals;
%        swarmFitnesses-(1*popSize). Fitness values of the individuals;
%Output: shelterIndexes-(1*dimension) vector. The updated shelter indexes
%                        of the individuals;
%        newShelterLeaders-(shelterNum*dimension) matrix. The updated
%                           shelter leaders;
%        newNormalizedLeaderFitnesses-(1*shelterNum) vector. The updated normalized
%                              fitness values of the shelter leaders;
%        swarmFitnesses-(1, popSize) vector. Updated fitness values of the
%                        individuals;
%Author: Fengji Luo
%Date: 07/2015
function [shelterIndexes, newShelterLeaders, newNormalizedLeaderFitnesses, swarmFitnesses] = updateShelters(dimension, swarm, shelterNum, shelterIndexes, swarmFitnesses)

[values, indexes] = sort(swarmFitnesses);  %firstly, sort the fitness values;


%second, set the shelter leaders;
newShelterLeaders = zeros(shelterNum, dimension); %the fitness function values of the imperialists;
for i=1:shelterNum
    originIndex = indexes(i);
    newShelterLeaders(i, :) = swarm(originIndex, :);
    shelterIndexes(1,originIndex) = i;
end

%third, get the leader fitness values;
newLeaderFitnesses = zeros(1, shelterNum); %the fitness function values of the imperialists;
for i=1:shelterNum
    newLeaderFitnesses(1, i) = values(i);
end

%fourth, set the normalized leader fitness values;
newNormalizedLeaderFitnesses = zeros(1, shelterNum);
for i=1:shelterNum
    newNormalizedLeaderFitnesses(1,i) = normalizeFitnessFunction(newLeaderFitnesses, newLeaderFitnesses(i), values(shelterNum+1));
end


end


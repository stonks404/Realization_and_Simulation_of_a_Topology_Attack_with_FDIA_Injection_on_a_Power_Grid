%Make the decision for an individual to leave its current shelter or not;
%Inputs: bounds-(2*dimension) matrix. The boundaries of each
%                dimension;
%        types-(1*dimension) vector. The data types of each dimension:
%        0-continous; 1-binary; 2-integer;
%        popSize-Integer. Number of the individuals;
%        dimension-Integer. Dimension of each individual;
%        shelterCaps-(1*shelterNum) vector. Capacities of the shelters;
%        shelterIndex-Integer. Index of the individual's current shelter;
%        theta-Real. Normalized fitness value of the shelter leader;
%        indNumberInShelter-Integer. Number of individuals present in the
%                           shelter;
%        shelterIndexes-(1*popSize) vector. Shelter indexes of the
%                        individuals;
%        swarm-(popSize*dimension) matrix. Set of the individuals;
%        indIndex-Integer. The index of the individual;
%        ShelterLeaders-(shelterNum*dimension) matrix. The set of the shelter leaders;
%        swarmFitnesses-(1*popSize). Fitness values of the individuals;
%        funcName_fitness-The name of the fitness function .m file;
%        userObj-Structure. Contain some extra information (if needed)
%                which are needed when calculate the fitness value; 
%        scale_local-Real. Scaling factor of the shelter leader;
%        alpha-Real. Movement factor;
%        Cr_local-Real. Local crossover factor;
%        Cr_global-Real. Global crossover factor;
%        bounceBack-Binary. 0-'let it fly'; 1-set the dimensional value to
%                   the boundary value when it exceeds the boundaries;
%Output: newSwarm-(popSize*dimension) matrix. The updated population;
%        newShelterIndexes-(1*popSize) vector. The updated shelter indexes
%                                              of the individuals;
%        newSwarmFitnesses-(1*popSize) vector. The updated fitness values
%                          of the individuals;
%Author: Fengji Luo
%Date: 07/2015
function [newSwarm, newShelterIndexes, newSwarmFitnesses] = decideLeave(bounds, types, popSize, dimension, shelterIndex, theta, indNumberInShelter, shelterIndexes, swarm, indIndex, shelterLeaders, swarmFitnesses,funcName_adjustInd, funcName_fitness, controlParams, userObj)   

shelterNum = controlParams.shelterNum;
shelterCaps = zeros(1, shelterNum) + controlParams.shelterCap;
scale_local = controlParams.scale_local;
Cr_local = controlParams.Cr_local;
alpha = controlParams.alpha;
Cr_global = controlParams.Cr_global;

bounceBack = controlParams.bounceBack;


%get the capacity of the given shelter;
shelterCap = shelterCaps(shelterIndex);

        
%calculate the probability for the cockroach to leave its corner to explore;
Q = calculateLeaveProbability(theta, shelterCap, indNumberInShelter);
leave = 0;
r = rand();
if(r<=Q)
    leave = 1;
    shelterIndexes(indIndex) = -1; %set the shelter index of the individual to be -1;
end

%do the search;
if(leave == 0)  %if the individual decides to stay in the shelter, do the local search;
    [newIndividual, newFitness] = localSearch(bounds, types, indIndex, dimension, swarm, shelterIndexes, shelterLeaders, swarmFitnesses, funcName_adjustInd, funcName_fitness, scale_local, Cr_local, bounceBack, userObj);
else %if the individual decides to leave the shelter, do the global search;
    [newIndividual, newFitness] = globalSearch(bounds, types, indIndex, dimension, popSize, swarm, swarmFitnesses, funcName_adjustInd, funcName_fitness, alpha, Cr_global, bounceBack, userObj);
end

swarm(indIndex, :) = newIndividual; %update the individual in the swarm;
swarmFitnesses(indIndex) = newFitness; %update the individual fitness;

%return the results;
newSwarm = swarm;
newShelterIndexes = shelterIndexes;
newSwarmFitnesses = swarmFitnesses;

end
%Make the decision for an individual to enter a given shelter or not;
%Inputs: bounds-(2*dimension) matrix. The boundaries of each
%                dimension;
%        types-(1*dimension) vector. The data types of each dimension:
%        0-continous; 1-binary; 2-integer;
%        popSize-Integer. Number of the individuals;
%        dimension-Integer. Dimension of each individual;
%        theta-Real. Normalized fitness value of the shelter leader;
%        shelterCaps-(1*shelterNum) vector. Capacities of the shelters;
%        shelterIndex-Integer. Index of the individual's current shelter;
%        indNumberInShelter-Integer. Number of individuals present in the
%                           shelter;
%        shelterNum-Integer. Number of the shelters;
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
function [newSwarm, newSubPopIndex, newFitnesses] = decideEnter(bounds, types, popSize, dimension, theta, shelterIndex, indNumberInShelter, shelterIndexes, swarm, indIndex, shelterLeader, swarmFitnesses, funcName_adjustInd, funcName_fitness, controlParams, userObj)

shelterNum = controlParams.shelterNum;
shelterCaps = zeros(1, shelterNum) + controlParams.shelterCap;
scale_local = controlParams.scale_local;
Cr_local = controlParams.Cr_local;
alpha = controlParams.alpha;
Cr_global = controlParams.Cr_global;

bounceBack = controlParams.bounceBack;


enter = 0;
if(shelterNum>0)
    S_corner = shelterCaps(shelterIndex);
    R = calculateEnterProbability(theta, S_corner, indNumberInShelter);
    r = rand();
    if(r<R)
        enter = 1;
        shelterIndexes(indIndex) = shelterIndex;  %decide to enter;
    end
end

         
if(enter == 0)
    [newCockroach, newFitness] = globalSearch(bounds, types, indIndex, dimension, popSize, swarm, swarmFitnesses, funcName_adjustInd, funcName_fitness, alpha, Cr_global, bounceBack, userObj);
    swarm(indIndex, :) = newCockroach;
    swarmFitnesses(indIndex) = newFitness;
else
    [newCockroach, newFitness] = localSearch(bounds, types, indIndex, dimension, swarm, shelterIndexes, shelterLeader, swarmFitnesses, funcName_adjustInd, funcName_fitness, scale_local, Cr_local, bounceBack, userObj);
    swarm(indIndex, :) = newCockroach;
    swarmFitnesses(indIndex) = newFitness;
end

newSwarm = swarm;
newSubPopIndex = shelterIndexes;
newFitnesses = swarmFitnesses;

end
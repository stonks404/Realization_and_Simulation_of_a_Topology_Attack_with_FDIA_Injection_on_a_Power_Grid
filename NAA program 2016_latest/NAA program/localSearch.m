%Do the local search;
%Inputs: bounds-(2*dimension) matrix. The boundaries of each
%                dimension;
%        types-(1*dimension) vector. The data types of each dimension:
%        0-continous; 1-binary; 2-integer;
%        indIndex-Integer. The index of the individual;
%        dimension-Integer. Dimension of each individual;
%        ShelterLeaders-(shelterNum*dimension) matrix. The set of the shelter leaders;
%        swarmFitnesses-(1*popSize). Fitness values of the individuals;
%        funcName_fitness-The name of the fitness function .m file;
%        userObj-Structure. Contain some extra information (if needed)
%                which are needed when calculate the fitness value; 
%        scale_local-Real. Scaling factor of the shelter leader;
%        Cr_local-Real. Local crossover factor;
%        bounceBack-Binary. 0-'let it fly'; 1-set the dimensional value to
%                   the boundary value when it exceeds the boundaries;
%Output: newInd-(1*dimension) vector. The generated individual after the local
%                              search;
%        newFitness-Real. The fitness value of the generated individual;
%Author: Fengji Luo
%Date: 07/2015
function [newInd, newFitness] = localSearch(bounds, types, indIndex, dimension, swarm, shelterIndexes, shelterLeaders, swarmFitnesses, funcName_adjustInd, funcName_fitness, scale_local, Cr_local, bounceBack, userObj)

individual = swarm(indIndex, :);
shelterIndex = shelterIndexes(indIndex);
indFitness = swarmFitnesses(indIndex);

newInd = individual;

shelterLeader = shelterLeaders(shelterIndex, :); %get the leader of current shelter;

%generate random number r;
int_r = unidrnd(dimension);

%crossover;
for i=1:dimension
     a = rand();
     if(individual == shelterLeader) %if it is the leader;
         if(a<Cr_local || i==int_r)
             lowerBound = -1*abs(individual(i)*scale_local);
             upperBound = abs(individual(i)*scale_local);
             newInd(i) = individual(i) + (lowerBound + rand()*(upperBound - lowerBound));
               
         end
     else  %if it is a follower;
        if(a<Cr_local || i==int_r)
            newInd(i) = individual(i)  + 2*rand() * (shelterLeader(i)-individual(i)); 
        end
     end
         
         
     %boundary constraint;
     if (bounceBack == 1)
         if (newInd(i)<bounds(1, i))
             newInd(i) = bounds(1, i) + (bounds(2,i)-bounds(1,i))*rand();
         elseif (newInd(i)>bounds(2, i))
             newInd(i) = bounds(1, i) + (bounds(2,i)-bounds(1,i))*rand();
         end
     end

end

%selection;
[newCockroach, newFitness] = evaluateCockroach(newInd, bounds, types, funcName_adjustInd, funcName_fitness, userObj);
newInd = newCockroach;
if(newFitness > indFitness)
   newInd = individual;
   newFitness = indFitness;
end


end
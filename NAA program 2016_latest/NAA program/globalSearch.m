%Do the local search;
%Inputs: bounds-(2*dimension) matrix. The boundaries of each
%                dimension;
%        types-(1*dimension) vector. The data types of each dimension:
%        0-continous; 1-binary; 2-integer;
%        indIndex-Integer. The index of the individual;
%        dimension-Integer. Dimension of each individual;
%        popSize-Integer. Number of the individuals;
%        swarm-(popSize*dimension) matrix. Set of the individuals;
%        swarmFitnesses-(1*popSize). Fitness values of the individuals;
%        funcName_fitness-The name of the fitness function .m file;
%        userObj-Structure. Contain some extra information (if needed)
%                which are needed when calculate the fitness value; 
%        alpha-Real. Movement factor;
%        Cr_global-Real. Global crossover factor;
%        bounceBack-Binary. 0-'let it fly'; 1-set the dimensional value to
%                   the boundary value when it exceeds the boundaries;
%Output: newInd-(1*dimension) vector. The generated individual after the global
%                              search;
%        newFitness-Real. The fitness value of the generated individual;
%Author: Fengji Luo
%Date: 07/2015
function [newInd, newFitness] = globalSearch(bounds, types, indIndex, dimension, popSize, swarm, swarmFitnesses, funcName_adjustInd, funcName_fitness, alpha, Cr_global, bounceBack, userObj)

individual = swarm(indIndex, :);
indFitness = swarmFitnesses(indIndex);
newInd = individual;

r_c_1 = unidrnd(popSize);
while(r_c_1==indIndex)
    r_c_1 = unidrnd(popSize);
end
r_c_2 = unidrnd(popSize);
while(r_c_2==indIndex || r_c_1 == r_c_2)
    r_c_2 = unidrnd(popSize);
end
r_c_3 = unidrnd(popSize);
while(r_c_3==indIndex || r_c_3 == r_c_2 || r_c_3 == r_c_1 )
    r_c_3 = unidrnd(popSize);
end

crossInd_1 = swarm(r_c_1, :);
crossInd_2 = swarm(r_c_2, :);
%generate random number r;
int_r = unidrnd(dimension);

%crossover;
for i=1:dimension
    a = rand();
    if( a<Cr_global || i==int_r)
        newInd(i) = individual(i) + alpha*rand() * (crossInd_1(i) - individual(i)) + alpha*rand()*(crossInd_2(i) - individual(i));
        
        %boundary constraint;
        if(bounceBack == 1)
             if (newInd(i)<bounds(1, i))
                 newInd(i) = bounds(1, i) + (bounds(2,i)-bounds(1,i))*rand();
             elseif (newInd(i)>bounds(2, i))
                 newInd(i) = bounds(1, i) + (bounds(2,i)-bounds(1,i))*rand();
             end
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




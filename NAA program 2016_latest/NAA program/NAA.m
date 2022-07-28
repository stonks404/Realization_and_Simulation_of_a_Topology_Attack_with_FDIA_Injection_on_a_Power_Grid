%Run NAA;
%Inputs: dimension-Integer. Dimension of the problem;
%        bounds-(2*dimension) matrix. The boundaries of each
%                dimension;
%        types-(1*dimension) matrix. The data type of each dimension:
%        0-continuous; 1-binary; 2-integer;
%        popSize-Integer. Number of the individuals;
%        iteration-Integer. Maximum iteration times;
%        funcName_adjustInd-The name of the individual adjustment function .m file;
%        funcName_fitness-The name of the fitness function .m file;
%        userObj-Structure. Contain some extra information (if needed)
%                which are needed when calculate the fitness value; 
%        scale_local-Real. Scaling factor of the shelter leader;
%                          controlParams-Structure. Control parameters of NAO. it includes
%                          following attributes: 
%                          [(1) shelterNum-Integer. number of the shelters;
%                           (2) shelterCap-Integer. capacity of each shelter; 
%                           (3) scale_local-Real. scaling factor of the shelter leaders; 
%                           (4) Cr_local-local crossover factor;
%                           (5) alpha-movement factor; 
%                           (6) Cr_global-global crossover factor;
%                           (7)bounceBack-boundary constraint handling flag. 0-'let it fly', 1-'bounce back'];
%        verbose-Iteration information display flag: 0-no display; 1-display;
%Outputs: bestFitness-Real. Found best fitness value; 
%         bestInd-(1*dimension) vector. The found best individual;
%         historicalBestFitness-(1*iteration) vector. The found best finess
%                                values at each iteration; 
%        historicalMeanFitness-(1*iteration) vector. The found mean finess
%                                values of the population at each iteration; 
%Author: Fengji Luo
%Date: 07/2015

function [bestFitness, bestInd, historicalBestFitness] = NAA(dimension, bounds, types, popSize, iteration, funcName_adjustInd, funcName_fitness, userObj, controlParams, verbose)

%control parameters;
shelterNum = controlParams.shelterNum;
shelterCap = zeros(1, shelterNum) + controlParams.shelterCap;


%initiallize the population;
[swarm] = initializeSwarm(popSize, dimension, bounds);  %initialize the individuals. Each individual represents a solution;
%initialize the shelters;
[swarm, swarmFitnesses, shelterIndexes, shelterLeaders, normalizedLeaderFitnesses]  = formShelters(swarm, dimension, bounds, types, popSize, shelterNum, shelterCap, funcName_adjustInd, funcName_fitness, userObj);

%initialize the global best fitness value and the corresponding individual;
[values_init, indexes_init] = sort(swarmFitnesses);
bestFitness = values_init(1);
bestInd = swarm(indexes_init(1), :);

% some statistically metrics;
historicalBestFitness = zeros(1, iteration);


%start the iteration process;
for i=1:iteration   
   indNumbers_shelters = calculateIndNumberInShelters(shelterNum, shelterIndexes);  %Calculate the number of individuals present in the shelters;
   
   %scan the individuals to do the search;
   for p=1:popSize
        shelterIndex = shelterIndexes(p);
        
        if(shelterIndex >0) %if the individual is currently in a shelter;
            theta = normalizedLeaderFitnesses(shelterIndex); %get the shelter quality value;
            [swarm, shelterIndexes, swarmFitnesses] = decideLeave(bounds, types, popSize, dimension, shelterIndex, theta, indNumbers_shelters(shelterIndex), shelterIndexes, swarm, p, shelterLeaders, swarmFitnesses, funcName_adjustInd, funcName_fitness, controlParams, userObj);
        else  %if the cockroach is currently doing the exploring;
            theta = 0;
            r_c = 0;
            numberInCorner = 0;
            if (shelterNum>0)
                r_c = unidrnd(shelterNum);
                theta = normalizedLeaderFitnesses(r_c);
                numberInCorner = indNumbers_shelters(r_c);
            end
            
            [swarm, shelterIndexes, swarmFitnesses] = decideEnter(bounds, types, popSize, dimension, theta, r_c, numberInCorner, shelterIndexes, swarm, p, shelterLeaders, swarmFitnesses, funcName_adjustInd, funcName_fitness, controlParams, userObj);
        end
        
   end

   
   %update the shelters;
   [shelterIndexes, shelterLeaders, normalizedLeaderFitnesses, swarmFitnesses] = updateShelters(dimension, swarm, shelterNum, shelterIndexes, swarmFitnesses);    
   
   
   %update the global best solution;
   [values_init, indexes_init] = sort(swarmFitnesses);
   if(values_init(1)<bestFitness)
       bestFitness = values_init(1);
       bestInd = swarm(indexes_init(1), :);
       indToEval = bestInd;
       for d=1:dimension
            value = bestInd(1, d);
            type = types(1, d);
            if(type==1)  %if the datatype of this dimension is binary;
                middle = (bounds(1, d) + bounds(2, d))/2;
                if(value<middle)
                    value = 0;
                else
                    value = 1;
                end
            elseif (type==2) %if the datatype of this dimension is integer;
                newValue = round(value);
                if(newValue < bounds(1, d))
                    newValue = ceil(bounds(1, d));
                elseif (newValue > bounds(2, d))
                    newValue = floor(bounds(2, d));
                end
       
                value = newValue;
            end
   
            indToEval(1, d) = value;
       end
       bestInd = indToEval;
   end
    historicalBestFitness(i) = bestFitness;
    
    
    
   if(verbose == 1)
        disp(['--------------------------------iteration ', num2str(i), ' best fitness:', num2str(bestFitness), '-----------------------------']);
   end

end

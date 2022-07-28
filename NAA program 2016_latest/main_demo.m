%
%the main file to invoke NAA;
%Author: Fengji Luo; Date: Aug/2016
%
%-----------------------------Prepare the parameters of NAA------------------------------------
%specify the population size and maximum generation time;
popSize = 40;
generation = 500;

%specify the dimension of the problem and bounds of each dimension;
D = 10;
bounds = [0, 0, 0, 0, 0, -10, -10, -10, 0, 0;
          100, 100, 100, 100, 100, 10, 10, 10, 1, 1];
      
      
%specify the data types of the problem: 0-continuous; 1-binary; 2-integer; 
types = [0,0,0,0,0,2,2,2,1,1];

%set the control parameters of NAA;
controlParam.shelterNum =2;
avg = popSize/(controlParam.shelterNum);
controlParam.shelterCap = avg;
controlParam.scale_local = 1;
controlParam.Cr_local =0.9; 
controlParam.Cr_global = 0.1;
controlParam.alpha = 1;

%specify whether to use the 'bounce back' strategy: 0-not use; 1-use;
%if choose to use the 'bounce back' strategy, then for a individual, if the
%value of any dimension exceeds the bounds, then a random number will be
%re-generated within the bounds;
%if choose not to use the 'bounce back' strategy, then no any adjustment
%will be applied on the individuals;
controlParam.bounceBack = 1;

%specify the file name of the fitness evaluation fuction;
fitnessFuncName = 'fitnessEval_demo';
adjustIndFuncName = 'constraintHandle_demo';

%specify the structure which contains the extra information used in the 
%constraint handling and fitness calculation modules;
userObj.threshold=80;

%------------------------------Run NAA-------------------------------------
% specify whether to output intermediate information:
%0-no display; 1-display
verbose = 1;

%run NAA.there are 3 outputs:
%bestFitness: the fitness value of the final found optimal solution;
%bestInd: the final found optimal solution;
%historicalFitness: a vector storing the fitness values of the best solution at each
%generation;

[bestFitness, bestInd, historicalFitness] = NAA(D, bounds, types, popSize,...
                                                generation, adjustIndFuncName, fitnessFuncName,...
                                                userObj, controlParam, verbose);
   

plot(historicalFitness);
xlabel('Generation');
ylabel('Fitness Value');
    












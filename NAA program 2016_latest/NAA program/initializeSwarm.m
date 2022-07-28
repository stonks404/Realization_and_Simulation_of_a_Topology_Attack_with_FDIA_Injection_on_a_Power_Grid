%initialize the population. Each individual represents a solution for the
%given problem.
%Inputs: popSize-Integer. The number of the individuals;
%        dimension-Integer. The dimension of each individual vector;
%        bounds-(2*dimension) matrix. entry(1, j) represents the initial
%               lower boundary of dimension j; entry (2, j) represents the
%               initial upper boundary of dimension j.
%Output: swarm-(popSize*dimension) matrix. Each row of it represents an
%               individual.
%Author: Fengji Luo
%Date: 07/2015
function [swarm] = initializeSwarm(popSize, dimension, bounds)

swarm = zeros(popSize, dimension); 

for i=1:popSize
    ind = zeros(1, dimension);
    for j=1:dimension
        ind(j) = bounds(1, j) + (bounds(2,j)-bounds(1,j))*rand();
    end 
    
    swarm(i, :) = ind;
end

end
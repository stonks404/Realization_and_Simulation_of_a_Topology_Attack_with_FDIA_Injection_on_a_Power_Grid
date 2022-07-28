%calculate the number of individuals present in each shelter;
%Inputs: shelterNum-Integer. Number of the shelters;
%        shelterIndexes-(1*popSize) vector. Shelter indexes of each
%                        individual;
%Output: numbers-(1*shelterNum) vector. Entry (1,i) represents the number
%                 of individuals present in the ith shelter;
%Author: Fengji Luo
%Date: 07/2015
function numbers = calculateIndNumberInShelters(shelterNum, shelterIndexes)

numbers = zeros(1, shelterNum);
popSize = length(shelterIndexes);

for p=1:popSize
    i = shelterIndexes(p);
    
    if(i>0)
        numbers(i) = numbers(i) + 1;
    end
end

end
%Calculate the probability of an expolore individual to enter a certain
%shelter;
%Inputs: theta-Real. Normalized fitness value of the candidate shelter
%                    leader;
%        shelterCap-Integer. Capacity of its current shelter;
%        indNumberInCorner-Integer. Number of individuals present in the
%                                   shelter;
%Output: R-Real. The entering probability;
%Author: Fengji Luo
%Date: 07/2015
function R = calculateEnterProbability(theta, shelterCap, indNumberInCorner)

R = (1-theta) * (1-(indNumberInCorner/shelterCap));

% R = 0.001 * (1-(indNumberInCorner/shelterCap));

end
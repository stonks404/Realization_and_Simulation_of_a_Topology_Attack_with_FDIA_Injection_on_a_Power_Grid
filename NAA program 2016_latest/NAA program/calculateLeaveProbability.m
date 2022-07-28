%Calculate the probability of an individual to leave its current shelter;
%Inputs: theta-Real. Normalized fitness value of its current shelter
%                    leader;
%        shelterCap-Integer. Capacity of its current shelter;
%        indNumberInCorner-Integer. Number of individuals present in the
%                                   shelter;
%Output: Q-Real. The leaving probability;
%Author: Fengji Luo
%Date: 07/2015
function Q = calculateLeaveProbability(theta, shelterCap, indNumberInCorner)

tmp = 1 + ( 1* ((indNumberInCorner/shelterCap)^2 ));
Q = theta / tmp;

%  Q = 0. / (1 + 1667*( (indNumberInCorner/shelterCap)^2 ) );


end
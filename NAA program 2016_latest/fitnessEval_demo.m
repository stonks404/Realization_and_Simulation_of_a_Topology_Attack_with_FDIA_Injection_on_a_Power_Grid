function fitness = fitnessEval_demo(ind, userObj)

dimension = length(ind);

fitness = 0;
for i=1:dimension
    v = ind(i);
    fitness = fitness + v*v;
end

end
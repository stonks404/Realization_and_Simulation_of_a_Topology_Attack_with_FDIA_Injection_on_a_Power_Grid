%Griewangk function;
function fitness = fitness_Griewangk(ind)

dimension = length(ind);

term1 = 0;
term2 = 1;
for i=1:dimension
    v = ind(i);
    term1 = term1 + (v*v/4000);
    term2 = term2 * cos(v/sqrt(i));
end
fitness = term1 - term2 + 1;

end
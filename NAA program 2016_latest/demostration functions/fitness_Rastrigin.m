%Rastrigin function;
function fitness = fitness_Rastrigin(ind)

dimension = length(ind);

fitness = dimension * 10;
for i=1:dimension
    v = ind(i);
    temp = v*v - 10 * cos(2*pi*v);
    fitness = fitness + temp;
end


end
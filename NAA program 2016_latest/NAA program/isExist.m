function [value, index] = isExist(array, element)
value = 0;
index = 0;

for i=1:length(array)
    if (array(1, i) == element)
        value = array(1,i);
        index = i;
        break;
    end
end

end
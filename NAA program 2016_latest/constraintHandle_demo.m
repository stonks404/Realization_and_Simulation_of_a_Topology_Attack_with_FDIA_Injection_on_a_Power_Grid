function [newInd, newUserObj] = constraintHandle_demo(ind, userObj)

limit = userObj.threshold;

if (ind(1)+ind(3)~=limit)
    if(ind(1)>limit)
        ind(1)=limit;
        ind(3)=0;
    else
        ind(3) = limit-ind(1);
    end
end

newInd = ind;
newUserObj = userObj;

end
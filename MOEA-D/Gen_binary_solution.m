function v = Gen_binary_solution(VarSize)
    cond = 0;
    v = randi([0 1], VarSize);
    while cond
        %display(cond);
        v = randi([0 1], VarSize);
        if (sum(v(1:33)) < 4)
            cond = 0;
        end
    end
end


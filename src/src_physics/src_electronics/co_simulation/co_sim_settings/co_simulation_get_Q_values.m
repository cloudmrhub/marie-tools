function[Q_val] = co_simulation_get_Q_values(RLC)

    Q_val = [];
    for i = 1:length(RLC)
        if RLC(i).optim.boolean == 1 || strcmp(RLC(i).type,'port')
            for j = 1:length(RLC(i).Q)
                Q_val = [Q_val RLC(i).Q(j)];
            end
        end
    end
    Q_val = Q_val.';

end
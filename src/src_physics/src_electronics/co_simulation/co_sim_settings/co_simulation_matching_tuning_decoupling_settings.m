function[all_matching_load_string,...
         all_port_order,...
         all_tuning_load_string,...
         all_tuning_mask_detailed,...
         all_matching_mask_detailed,...
         all_fixed_matching_mask_detailed,...
         all_fixed_matching_values,...
         all_matching_mask,...
         all_fixed_matching_mask,...
         all_tuning_mask,...
         all_mutual_mask] = co_simulation_matching_tuning_decoupling_settings(matching_load_string,...
                                                                              tuning_load_string,...
                                                                              port_order,...
                                                                              tuning_mask_detailed,...
                                                                              mutual_mask_detailed,...
                                                                              matching_mask_detailed,...
                                                                              fixed_matching_mask_detailed,...
                                                                              fixed_matching_values,...
                                                                              matching_mask,...
                                                                              fixed_matching_mask,...
                                                                              tuning_mask,...
                                                                              mutual_mask)
   
    max_len = max(cellfun(@(x) numel(x), {port_order.id}));
    n_ports = numel(port_order);
    all_port_order = zeros(max_len, n_ports);
    for i = 1:n_ports
        ids = port_order(i).id;
        all_port_order(1:numel(ids), i) = ids;
    end
    all_matching_load_string = strings(size(all_port_order));  
    for i = 1:numel(matching_load_string)
        ids = matching_load_string(i).id;
        all_matching_load_string(1:numel(ids), i) = ids;
    end
    
    mask = (tuning_mask_detailed+mutual_mask_detailed) ~= 0;
    [N, M] = size(tuning_mask_detailed);
    str_matrix = strings(N, M);  
    for i = 1:M
        str_matrix(mask(:,i),i) = tuning_load_string(i).id;
    end
    all_tuning_load_string = strings(N, 1);
    for i = 1:N
        row_vals = str_matrix(i, :);
        nonempty = row_vals ~= "";
        if any(nonempty)
            all_tuning_load_string(i) = row_vals(find(nonempty, 1, 'first'));
        end
    end
    all_tuning_load_string(all_tuning_load_string == "") = [];


    all_fixed_matching_values                                                 = sum(fixed_matching_values, 2) ./ sum(fixed_matching_values ~= 0, 2);
    all_fixed_matching_values(isnan(all_fixed_matching_values))               = 0;
    all_matching_mask_detailed                                                = sum(matching_mask_detailed, 2) ./ sum(matching_mask_detailed ~= 0, 2);
    all_matching_mask_detailed(isnan(all_matching_mask_detailed))             = 0;
    all_fixed_matching_mask_detailed                                          = sum(fixed_matching_mask_detailed, 2) ./ sum(fixed_matching_mask_detailed ~= 0, 2);
    all_fixed_matching_mask_detailed(isnan(all_fixed_matching_mask_detailed)) = 0;
    all_tuning_mask_detailed                                                  = sum(tuning_mask_detailed, 2) ./ sum(tuning_mask_detailed ~= 0, 2);
    all_tuning_mask_detailed(isnan(all_tuning_mask_detailed))                 = 0;
    all_matching_mask                                                         = sum(matching_mask, 2) ./ sum(matching_mask ~= 0, 2);
    all_matching_mask(isnan(all_matching_mask))                               = 0;
    all_fixed_matching_mask                                                   = sum(fixed_matching_mask, 2) ./ sum(fixed_matching_mask ~= 0, 2);
    all_fixed_matching_mask(isnan(all_fixed_matching_mask))                   = 0;
    all_tuning_mask                                                           = sum(tuning_mask, 2) ./ sum(tuning_mask ~= 0, 2);
    all_tuning_mask(isnan(all_tuning_mask))                                   = 0;
    all_mutual_mask                                                           = sum(mutual_mask, 3) ./ sum(mutual_mask ~= 0, 3);
    all_mutual_mask(isnan(all_mutual_mask))                                   = 0;

end
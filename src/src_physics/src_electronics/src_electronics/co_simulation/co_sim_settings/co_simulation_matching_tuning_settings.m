function[symmetries_all,...
         new_port_order,...
         new_loads_matching,...
         new_loads_tuning,...
         new_symms,...
         matching_mask_detailed,...
         fixed_matching_mask_detailed,...
         fixed_matching_values,...
         tuning_mask_detailed,...
         mutual_mask_detailed,...
         lower_bound_all,...
         upper_bound_all,...
         unique_indices_all,...
         unique_indices] = co_simulation_matching_tuning_settings(RLC_tmd,...
                                                                  RLC_org,...
                                                                  tuning_elements,...
                                                                  mutual_elements,...
                                                                  elems_after_tuning,...
                                                                  coil_entity,...
                                                                  fval_tune)
    

    max_entities = max(coil_entity,[],'all');

    matching_element_counter = 0;
    fixed_matching_element_counter = 0;
    for i = 1:length(RLC_tmd)
        if RLC_tmd(i).optim.boolean == 1 && strcmp(RLC_org(i).type,'port')
            matching_element_counter = matching_element_counter + length(convertCharsToStrings(RLC_tmd(i).load));
        elseif RLC_tmd(i).optim.boolean == 0 && strcmp(RLC_org(i).type,'port')
            fixed_matching_element_counter = fixed_matching_element_counter + length(convertCharsToStrings(RLC_tmd(i).load));
        end
    end
    all_matching_element_counter = matching_element_counter + fixed_matching_element_counter;

    % Set up optimization arrays
    loads_matching               = strings(all_matching_element_counter,max_entities);
    port_order                   = zeros(all_matching_element_counter,max_entities);
    fixed_matching_values        = zeros(all_matching_element_counter+tuning_elements+mutual_elements,max_entities);
    loads_tuning                 = strings(tuning_elements+mutual_elements,max_entities);
    symmetries_all               = zeros(all_matching_element_counter+tuning_elements+mutual_elements,1);
    symmetries                   = zeros(matching_element_counter+tuning_elements+mutual_elements,max_entities);
    lower_bound_all              = nan(all_matching_element_counter+tuning_elements+mutual_elements,1);
    upper_bound_all              = nan(all_matching_element_counter+tuning_elements+mutual_elements,1);
    matching_mask_detailed       = zeros(all_matching_element_counter+tuning_elements+mutual_elements,max_entities);
    fixed_matching_mask_detailed = zeros(all_matching_element_counter+tuning_elements+mutual_elements,max_entities);
    tuning_mask_detailed         = zeros(all_matching_element_counter+tuning_elements+mutual_elements,max_entities);
    mutual_mask_detailed         = zeros(all_matching_element_counter+tuning_elements+mutual_elements,max_entities);

    c1 = 0;
    c2 = 0;
    c3 = 0;
    t  = 0;
    for port_num = 1:length(RLC_org)
        port_order_counter = 1;
        if RLC_tmd(port_num).optim.boolean == 1 || (RLC_tmd(port_num).optim.boolean == 0 && strcmp(RLC_org(port_num).type,'port'))
            if strcmp(RLC_org(port_num).type,'port') && RLC_tmd(port_num).optim.boolean == 1 
                for i = 1:length(convertCharsToStrings(RLC_org(port_num).load))
                    load_tmp = convertCharsToStrings(RLC_tmd(port_num).load);
                    for j = 1:size(coil_entity,2)
                        if coil_entity(port_num,j)>0
                            matching_mask_detailed(c1+1,coil_entity(port_num,j)) = c1+1;
                            loads_matching(c3+1,coil_entity(port_num,j))         = load_tmp(i);
                            port_order(c3+1,coil_entity(port_num,j))             = port_order_counter;
                            symmetries(c1+1,coil_entity(port_num,j))             = RLC_tmd(port_num).optim.symmetry+c1;
                        end
                    end
                    symmetries_all(c1+1)  = RLC_tmd(port_num).optim.symmetry+c1;
                    lower_bound_all(c1+1) = RLC_tmd(port_num).optim.minim(i);
                    upper_bound_all(c1+1) = RLC_tmd(port_num).optim.maxim(i);
                    c1 = c1+1;
                    c3 = c3+1;
                    port_order_counter = port_order_counter+1;
                end
                t  = c1-1;
            elseif strcmp(RLC_org(port_num).type,'port') && RLC_tmd(port_num).optim.boolean == 0
                for i = 1:length(convertCharsToStrings(RLC_org(port_num).load))
                    load_tmp = convertCharsToStrings(RLC_tmd(port_num).load);
                    value_tmp = RLC_tmd(port_num).value(i);
                    for j = 1:size(coil_entity,2)
                        if coil_entity(port_num,j)>0
                            fixed_matching_mask_detailed(c1+1,coil_entity(port_num,j)) = c1+1;
                            loads_matching(c3+1,coil_entity(port_num,j))               = load_tmp(i);
                            fixed_matching_values(c1+1,coil_entity(port_num,j))        = value_tmp;
                            port_order(c3+1,coil_entity(port_num,j))                   = port_order_counter;
                        end
                    end
                    symmetries_all(c1+1) = RLC_tmd(port_num).optim.symmetry+c1;
                    c1 = c1+1;
                    c3 = c3+1;
                    port_order_counter = port_order_counter+1;
                end
                t  = c1-1;
            elseif strcmp(RLC_org(port_num).type,'element')
                c2 = c2 + 1;
                for j = 1:size(coil_entity,2)
                    if coil_entity(port_num,j)>0
                        tuning_mask_detailed(c1+1,coil_entity(port_num,j)) = c1+1;
                        loads_tuning(c2,coil_entity(port_num,j))           = convertCharsToStrings(RLC_tmd(port_num).load);
                        symmetries(c1+1,coil_entity(port_num,j))           = RLC_tmd(port_num).optim.symmetry+t;
                    end
                end
                symmetries_all(c1+1)  = RLC_tmd(port_num).optim.symmetry+t;

                if fval_tune(coil_entity(port_num,1)) > 1
                    lower_bound_all(c1+1) = RLC_tmd(port_num).optim.minim;
                    upper_bound_all(c1+1) = RLC_tmd(port_num).optim.maxim;
                else
                    lower_bound_all(c1+1) = 0.8*elems_after_tuning(c2);
                    upper_bound_all(c1+1) = 1.2*elems_after_tuning(c2);
                end
                c1 = c1+1;
            end
        end
    end    

    mutual_mod = 1;
    large_symmetry_counter = 1;
    for port_num = 1:length(RLC_org)
        if RLC_tmd(port_num).optim.boolean == 1 && strcmp(RLC_org(port_num).type,'element') && strcmp(RLC_org(port_num).load,'mutual_inductor')
            if mod(mutual_mod,2)
                c2 = c2 + 1;
                for j = 1:size(coil_entity,2)
                    if coil_entity(port_num,j)>0
                        loads_tuning(c1+1,coil_entity(port_num,j))         = convertCharsToStrings('mutual_inductance_coefficient');
                        mutual_mask_detailed(c1+1,coil_entity(port_num,j)) = c1+1;
                        symmetries(c1+1,coil_entity(port_num,j))           = large_symmetry_counter+matching_element_counter+tuning_elements+mutual_elements;
                    end
                end
                symmetries_all(c1+1) = large_symmetry_counter+matching_element_counter+tuning_elements+mutual_elements;
                if elems_after_tuning(c2) > 0
                    lower_bound_all(c1+1) = 0.8*elems_after_tuning(c2);
                    upper_bound_all(c1+1) = min(1,1.2*elems_after_tuning(c2));
                else
                    lower_bound_all(c1+1) = max(-1,1.2*elems_after_tuning(c2));
                    upper_bound_all(c1+1) = 0.8*elems_after_tuning(c2);
                end
                large_symmetry_counter    = large_symmetry_counter + 1;
                c1 = c1+1;
            end
            mutual_mod = mutual_mod + 1;
        end
    end  

    [unique_vals, unique_indices_all, new_indices] = unique(nonzeros(symmetries_all), 'stable'); 
    new_mapping                                    = 1:length(unique_vals);
    symmetries_all                                 = new_mapping(new_indices);

    unique_indices     = repmat(struct('id',[]), max_entities, 1);
    new_symms          = repmat(struct('id',[]), max_entities, 1);
    new_loads_matching = repmat(struct('id',[]), max_entities, 1);
    new_loads_tuning   = repmat(struct('id',[]), max_entities, 1);
    new_port_order     = repmat(struct('id',[]), max_entities, 1);
    for i = 1:max_entities
        q                                        = symmetries(:,i);
        nonzero_idx                              = find(q ~= 0);
        nonzero_vals                             = q(nonzero_idx);
        [unique_vals, ia_nonzero, new_symms_tmp] = unique(nonzero_vals, 'stable');
        new_mapping                              = 1:length(unique_vals);
        unique_indices(i).id                     = nonzero_idx(ia_nonzero);
        new_symms(i).id                          = new_mapping(new_symms_tmp);
        new_loads_matching(i).id                 = loads_matching(strlength(loads_matching(:,i))>0,i);
        new_loads_tuning(i).id                   = loads_tuning(strlength(loads_tuning(:,i))>0,i);
        new_port_order(i).id                     = nonzeros(port_order(:,i));
    end

end
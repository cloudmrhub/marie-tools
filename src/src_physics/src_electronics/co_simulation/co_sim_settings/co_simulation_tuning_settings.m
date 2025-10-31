function[tuning_elements,...
         mutual_elements,...
         new_loads,...
         new_symms,...
         matching_mask,...
         fixed_matching_mask,...
         tuning_mask,...
         mutual_mask,...
         unique_indices,...
         lower_bound,...
         upper_bound,...
         coil_entity] = co_simulation_tuning_settings(RLC_tmd,RLC_org)

    % Get variable elements
    tuning_elements         = 0;
    matching_elements       = 0;
    mutual_elements         = 0;
    fixed_matching_elements = 0;
    for port_num = 1:length(RLC_org)
        if strcmp(RLC_org(port_num).type,'port')
            if RLC_tmd(port_num).optim.boolean == 1
                matching_elements = matching_elements+1;
            else
                fixed_matching_elements = fixed_matching_elements+1;
            end
        elseif RLC_tmd(port_num).optim.boolean == 1 && strcmp(RLC_org(port_num).type,'element')
            tuning_elements = tuning_elements+1;
            if strcmp(RLC_org(port_num).load,'mutual_inductor')
                mutual_elements = mutual_elements+1/2;
            end
        end
    end
    mutual_elements = round(mutual_elements);
    all_matching_elements = fixed_matching_elements+matching_elements;
    
    all_elements = all_matching_elements+tuning_elements+mutual_elements;

    coil_entity = zeros(all_matching_elements+tuning_elements,all_matching_elements);
    counter_ports = 1;
    for port_num = 1:length(RLC_org)
        if RLC_tmd(port_num).optim.boolean == 1 || (strcmp(RLC_org(port_num).type,'port') && RLC_tmd(port_num).optim.boolean == 0)
            coil_entity(counter_ports,1:length(RLC_org(port_num).excitation.entity)) = RLC_org(port_num).excitation.entity;
            counter_ports = counter_ports + 1;
        end
    end
    max_entities = max(coil_entity,[],'all');
    coil_entity = coil_entity(:, any(coil_entity, 1));

    % Set up optimization arrays
    loads_all           = strings(all_elements,max_entities);
    symmetries_tun      = zeros(tuning_elements+mutual_elements,max_entities);
    lower_bound         = zeros(all_elements,1);
    upper_bound         = zeros(all_elements,1);
    matching_mask       = zeros(all_elements,max_entities);
    fixed_matching_mask = zeros(all_elements,max_entities);
    tuning_mask         = zeros(all_elements,max_entities);
    mutual_mask         = zeros(all_elements,2,max_entities);
    c1 = 0;
    c2 = 0;
    c3 = 0;
    c4 = 0;
    c5 = 0;
    mutual_mod = 1;
    for port_num = 1:length(RLC_org)
        if RLC_tmd(port_num).optim.boolean == 1 || (strcmp(RLC_org(port_num).type,'port') && RLC_tmd(port_num).optim.boolean == 0)
            c1 = c1+1;
            if strcmp(RLC_org(port_num).type,'port') && RLC_tmd(port_num).optim.boolean == 1
                c2 = c2+1;
                for j = 1:size(coil_entity,2)
                    if coil_entity(port_num,j)>0
                        matching_mask(c2,coil_entity(port_num,j)) = c1;
                    end
                end
            elseif strcmp(RLC_org(port_num).type,'port') && RLC_tmd(port_num).optim.boolean == 0
                c5 = c5+1;
                for j = 1:size(coil_entity,2)
                    if coil_entity(port_num,j)>0
                        fixed_matching_mask(c5,coil_entity(port_num,j)) = c1;
                    end
                end
            elseif strcmp(RLC_org(port_num).type,'element')
                c3 = c3+1;
                for j = 1:size(coil_entity,2)
                    if coil_entity(port_num,j)>0
                        tuning_mask(c3,coil_entity(port_num,j)) = c1;
                    end
                end
                if strcmp(RLC_org(port_num).load,'mutual_inductor')
                    if mod(mutual_mod,2)
                        c4 = c4+1;
                        for j = 1:size(coil_entity,2)
                            if coil_entity(port_num,j)>0
                                mutual_mask(c4,1,coil_entity(port_num,j)) = RLC_org(RLC_org(port_num).cross_talk(1)).cross_talk(1);
                                mutual_mask(c4,2,coil_entity(port_num,j)) = RLC_org(port_num).cross_talk(1);
                            end
                        end
                    end
                    mutual_mod = mutual_mod + 1;
                end
            end
        end
    end 

    c1 = 0;
    for port_num = 1:length(RLC_org)
        if RLC_tmd(port_num).optim.boolean == 1 && strcmp(RLC_org(port_num).type,'element')
            c1 = c1+1;
            for j = 1:size(coil_entity,2)
                if coil_entity(port_num,j)>0
                    loads_all(c1,coil_entity(port_num,j))      = convertCharsToStrings(RLC_tmd(port_num).load);
                    symmetries_tun(c1,coil_entity(port_num,j)) = RLC_tmd(port_num).optim.symmetry;
                end
            end
            lower_bound(c1,1) = RLC_tmd(port_num).optim.minim;
            upper_bound(c1,1) = RLC_tmd(port_num).optim.maxim;
        end
    end    

    mutual_mod = 1;
    large_symmetry_counter = 1;
    for port_num = 1:length(RLC_org)
        if RLC_tmd(port_num).optim.boolean == 1 && strcmp(RLC_org(port_num).type,'element') && strcmp(RLC_org(port_num).load,'mutual_inductor')
            if mod(mutual_mod,2)
                c1 = c1+1;
                for j = 1:size(coil_entity,2)
                    if coil_entity(port_num,j)>0
                        loads_all(c1,coil_entity(port_num,j))      = convertCharsToStrings('mutual_inductance_coefficient');
                        symmetries_tun(c1,coil_entity(port_num,j)) = large_symmetry_counter+all_elements;
                    end
                end
                k1                     = RLC_org(port_num).cross_talk(2)/sqrt(RLC_tmd(port_num).optim.maxim*RLC_tmd(RLC_org(port_num).cross_talk(1)).optim.maxim);
                k2                     = RLC_org(port_num).cross_talk(2)/sqrt(RLC_tmd(port_num).optim.minim*RLC_tmd(RLC_org(port_num).cross_talk(1)).optim.minim);
                lower_bound(c1,1)      = min(max(min(k1, 1), -1), max(min(k2, 1), -1));
                upper_bound(c1,1)      = max(max(min(k1, 1), -1), max(min(k2, 1), -1));
                large_symmetry_counter = large_symmetry_counter + 1;
            end
            mutual_mod = mutual_mod + 1;
        end
    end    

    unique_indices = repmat(struct('id',[]), max_entities, 1);
    new_symms      = repmat(struct('id',[]), max_entities, 1);
    new_loads      = repmat(struct('id',[]), max_entities, 1);
    for i = 1:max_entities
        q                                        = symmetries_tun(:,i);
        nonzero_idx                              = find(q ~= 0);
        nonzero_vals                             = q(nonzero_idx);
        [unique_vals, ia_nonzero, new_symms_tmp] = unique(nonzero_vals, 'stable');
        new_mapping                              = 1:length(unique_vals);
        unique_indices(i).id                     = nonzero_idx(ia_nonzero);
        new_symms(i).id                          = new_mapping(new_symms_tmp);
        new_loads(i).id = loads_all(strlength(loads_all(:,i))>0,i);
    end

end
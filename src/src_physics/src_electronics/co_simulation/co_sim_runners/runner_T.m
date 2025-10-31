function[elems_after_tuning,...
         tuning_elements,...
         mutual_elements,...
         tuning_mask,...
         matching_mask,...
         fixed_matching_mask,...
         mutual_mask,...
         coil_entity,...
         fval_min_return] = runner_T(RLC_tmd,...
                                     RLC_org,...
                                     options_PSO_1,...
                                     YPn,...
                                     emc)

    [tuning_elements,...
    mutual_elements,...
    tuning_load_string,...
    new_symms,...
    matching_mask,...
    fixed_matching_mask,...
    tuning_mask,...
    mutual_mask,...
    unique_indices,...
    lower_bound,...
    upper_bound,...
    coil_entity] = co_simulation_tuning_settings(RLC_tmd,RLC_org);
    
    elems_after_tuning = zeros(length(lower_bound),1);
    fval_min_return = inf*ones(max(coil_entity,[],'all'),1);

    if ~all(tuning_mask(:) == 0)

        for i = 1:max(coil_entity,[],'all')

            if isempty(unique_indices(i).id)
                continue
            end
        
            [y_ids,~]       = find(coil_entity == i);
            non_ids         = setdiff(1:size(YPn, 1), y_ids);
            YPn_to_pass     = YPn(y_ids,y_ids)-YPn(y_ids,non_ids)/YPn(non_ids,non_ids)*YPn(non_ids,y_ids);
            M               = nonzeros(matching_mask(:,i));
            F               = nonzeros(fixed_matching_mask(:,i));
            N               = nonzeros(tuning_mask(:,i));
            D1              = nonzeros(mutual_mask(:,1,i));
            D2              = nonzeros(mutual_mask(:,2,i));
            D1              = ismember(N,D1);
            D2              = ismember(N,D2);
            [M,F,N]         = map_to_local(M,F,N);
            FM              = sort([M(:); F(:)]);
            D1              = N(D1);
            D2              = N(D2);
            all_values      = [D1 D2];                     
            [~, sort_idx]   = sort(all_values);  
            ranks           = zeros(size(all_values));           
            ranks(sort_idx) = 1:length(all_values); 
            mapped_D1       = ranks(1:length(D1));
            mapped_D2       = ranks(length(D1)+1:end);
            idx12           = sub2ind(size(YPn_to_pass), D1, D2);
            idx21           = sub2ind(size(YPn_to_pass), D2, D1);

            costFun_tune = @(X) tuning(X,...
                                       YPn_to_pass,...
                                       new_symms(i).id,...
                                       FM,...
                                       N,...
                                       idx12,...
                                       idx21,...
                                       mapped_D1,...
                                       mapped_D2,...
                                       tuning_load_string(i).id,...
                                       emc.omega);
            
            
            fval_min = inf;
            for counter_repeat = 1:2

                rng(counter_repeat);
                
                [X_final, fval, ~, ~] = particleswarm(costFun_tune,...
                                                      length(unique_indices(i).id),...
                                                      lower_bound(unique_indices(i).id),...
                                                      upper_bound(unique_indices(i).id),...
                                                      options_PSO_1);
                if fval<fval_min
                    fval_min = fval;
                    X_out = X_final;
                    options_PSO_1.SwarmSize = min(options_PSO_1.SwarmSize+100,1000);
                    options_PSO_1.MaxIterations = min(options_PSO_1.MaxIterations+100,1500);
                    options_PSO_1.InertiaRange = [max(0.3,options_PSO_1.InertiaRange(1)-0.1) min(options_PSO_1.InertiaRange(2)+0.1,2.0)];
                    options_PSO_1.SelfAdjustmentWeight = max(0.2,options_PSO_1.SelfAdjustmentWeight-0.1);
                    options_PSO_1.SocialAdjustmentWeight = min(options_PSO_1.SocialAdjustmentWeight+0.1,2.5);
                end
    
                if fval_min < 1
                    break
                end
    
            end
            
            if isempty(unique_indices(i).id)
                break
            end

            if fval_min > 1
                fprintf('\t Coil %i is not tuned, %4.6f.\n',i,fval_min);
            else
                fprintf('\t Coil %i is tuned, %4.6f.\n',i,fval_min);
            end
    
            % Update in case elements are shared
            Xmin_raw                          = 0.5 * X_out;
            Xmax_raw                          = 1.5 * X_out;
            X_min_fixed                       = max(min(min(Xmin_raw, Xmax_raw), 1), -1);
            X_max_fixed                       = max(min(max(Xmin_raw, Xmax_raw), 1), -1);
            lower_bound(unique_indices(i).id) = X_min_fixed;
            upper_bound(unique_indices(i).id) = X_max_fixed;
    
            id_tune = nonzeros(tuning_mask(:,i));
            id_mutual = find(mutual_mask(:, 1, i) ~= 0) + max(tuning_mask,[],'all');
            id_store = [id_tune; id_mutual];
            elems_after_tuning(id_store) = X_out(new_symms(i).id);
    
        end
        fval_min_return(i) = fval_min;
        elems_after_tuning = nonzeros(elems_after_tuning);

    else

        elems_after_tuning = [];
    
    end

end
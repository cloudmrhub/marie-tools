function[X_out,...
         all_matching_load_string,...
         all_port_order,...
         all_tuning_load_string,...
         all_tuning_mask_detailed,...
         all_matching_mask_detailed,...
         all_fixed_matching_mask_detailed,...
         all_fixed_matching_values,...
         all_matching_mask,...
         all_fixed_matching_mask,...
         all_tuning_mask,...
         all_mutual_mask,...
         match_split,...
         mask_detailed_split,...
         fval_min,...
         check_Rx,...
         check_Tx,...
         check_TxRx] = runner_M_T_D_PD_split(RLC_tmd,...
                                             RLC_org,...
                                             matching_load_string,...
                                             tuning_load_string,...
                                             port_order,...
                                             symmetries_all,...
                                             matching_mask_detailed,...
                                             fixed_matching_mask_detailed,...
                                             fixed_matching_values,...
                                             tuning_mask_detailed,...
                                             mutual_mask_detailed,...
                                             matching_mask,...
                                             fixed_matching_mask,...
                                             tuning_mask,...
                                             mutual_mask,...
                                             unique_indices_all,...
                                             lower_bound_mat,...
                                             upper_bound_mat,...
                                             elems_after_matching,...
                                             options_PSO_2,...
                                             YPn,...
                                             emc,...
                                             fval_match_all)

    [all_matching_load_string,...
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
                                                                          mutual_mask);

    [match_split,...
     tune_split,...
     mask_split,...
     matching_mask_split,...
     fixed_matching_mask_split,...
     tuning_mask_split,...
     mutual_mask_split,...
     mask_detailed_split,...
     fixed_mask_detailed_split] = co_simulation_split_settings(RLC_tmd,...
                                                               RLC_org);

    Rx_flag = ~isempty(all_matching_load_string(:,match_split==1));
    Tx_flag = ~isempty(all_matching_load_string(:,match_split==2));
    TxRx_flag = ~isempty(all_matching_load_string(:,match_split==3));

    Rx_elems_detailed   = mask_detailed_split==1;
    Tx_elems_detailed   = mask_detailed_split==2;
    TxRx_elems_detailed = mask_detailed_split==3;
    
    if Rx_flag
        Rx_elems_match                                          = match_split==1;
        Rx_elems_tune                                           = tune_split==1;
        Rx_elems_matching_mask_split                            = matching_mask_split==1;
        Rx_elems_fixed_matching_mask_split                      = fixed_matching_mask_split==1;
        Rx_elems_tuning_mask_split                              = tuning_mask_split==1;
        Rx_elems_mutual_mask_split                              = mutual_mask_split==1;
        symmetries_all_to_pass                                  = symmetries_all(Rx_elems_detailed)';
        [~, ~, symmetries_all_to_pass]                          = unique(symmetries_all_to_pass, 'stable'); 
        
        fval_match = sum(fval_match_all(Rx_elems_match));
        if fval_match > 1
            low_lv = 0.5;
            upp_lv = 1.5;
        else
            low_lv = 0.8;
            upp_lv = 1.2;
        end
        if ~isempty(nonzeros(Rx_elems_detailed))

            [y_ids,~]   = find(mask_split == 1);
            non_ids     = setdiff(1:size(YPn, 1), y_ids);
            YPn_to_pass = YPn(y_ids,y_ids)-YPn(y_ids,non_ids)/YPn(non_ids,non_ids)*YPn(non_ids,y_ids);

            M                                  = nonzeros(all_matching_mask(Rx_elems_matching_mask_split));
            F                                  = nonzeros(all_fixed_matching_mask(Rx_elems_fixed_matching_mask_split));
            N                                  = nonzeros(all_tuning_mask(Rx_elems_tuning_mask_split));
            D1                                 = nonzeros(all_mutual_mask(Rx_elems_mutual_mask_split,1));
            D2                                 = nonzeros(all_mutual_mask(Rx_elems_mutual_mask_split,2));
            D1                                 = ismember(N,D1);
            D2                                 = ismember(N,D2);
            [M,F,N]                            = map_to_local(M,F,N);
            FM                                 = sort([M(:); F(:)]);
            [M_new,~]                          = map_to_local_2(M,F);
            D1                                 = N(D1);
            D2                                 = N(D2);
            all_values                         = [D1 D2];                     
            [~, sort_idx]                      = sort(all_values);  
            ranks                              = zeros(size(all_values));           
            ranks(sort_idx)                    = 1:length(all_values); 
            mapped_D1                          = ranks(1:length(D1));
            mapped_D2                          = ranks(length(D1)+1:end);
            length_matching_network            = unique(max(all_port_order(:,Rx_elems_match)));
            idx12                              = sub2ind(size(YPn), D1, D2);
            idx21                              = sub2ind(size(YPn), D2, D1);
            M_detailed                         = nonzeros(all_matching_mask_detailed(Rx_elems_detailed));
            F_detailed                         = nonzeros(all_fixed_matching_mask_detailed(Rx_elems_detailed));
            N_detailed                         = nonzeros(all_tuning_mask_detailed(Rx_elems_detailed));
            [M_detailed,F_detailed,N_detailed] = map_to_local(M_detailed,F_detailed,N_detailed);
            MN_detailed                        = sort([M_detailed(:); N_detailed(:)]);
            MF_detailed                        = sort([M_detailed(:); F_detailed(:)]);
            symmetries_all_to_pass(F_detailed) = [];
            [unique_vals, ~, new_indices]      = unique(nonzeros(symmetries_all_to_pass), 'stable'); 
            new_mapping                        = 1:length(unique_vals);
            symmetries_all_to_pass             = new_mapping(new_indices);

            costFun_match_tune_preamplifierdecoupling = @(X) matching_and_tuning_and_preamplifierdecoupling(X,...
                                                                                                            YPn_to_pass,...
                                                                                                            symmetries_all_to_pass,...
                                                                                                            M_detailed,...
                                                                                                            F_detailed,...
                                                                                                            N_detailed,...
                                                                                                            MN_detailed,...
                                                                                                            MF_detailed,...
                                                                                                            all_fixed_matching_values(Rx_elems_detailed),...
                                                                                                            FM,...
                                                                                                            M_new,...
                                                                                                            N,...
                                                                                                            idx12,...
                                                                                                            idx21,...
                                                                                                            mapped_D1,...
                                                                                                            mapped_D2,...
                                                                                                            length_matching_network,...
                                                                                                            all_matching_load_string(:,Rx_elems_match),...
                                                                                                            all_tuning_load_string(Rx_elems_tune),...
                                                                                                            all_port_order(:,Rx_elems_match),...
                                                                                                            emc.omega,...
                                                                                                            emc.Z0,...
                                                                                                            emc.Preamp_res);
            fval_min_Rx = inf;
            if fval_match > 50
                lb_t = lower_bound_mat(unique_indices_all);
                ub_t = upper_bound_mat(unique_indices_all);
                lub_m = Rx_elems_detailed(unique_indices_all);
                X_min_fixed = rmmissing(lb_t(lub_m));
                X_max_fixed = rmmissing(ub_t(lub_m));
            else
                lub_m = Rx_elems_detailed(unique_indices_all);
                Xmin_raw = low_lv * rmmissing(elems_after_matching(lub_m));
                Xmax_raw = upp_lv * rmmissing(elems_after_matching(lub_m));
                X_min_fixed = max(min(min(Xmin_raw, Xmax_raw), 1), -1);
                X_max_fixed = max(min(max(Xmin_raw, Xmax_raw), 1), -1);
            end
    
            for counter_repeat = 1:10
                if isempty(X_min_fixed)
                    X_out_Rx = [];
                    break;
                end
                if fval_min_Rx < 1
                    break
                end
                rng(counter_repeat);
    
                [X_final, fval, ~, ~] = particleswarm(costFun_match_tune_preamplifierdecoupling,...
                                                      length(X_min_fixed),...
                                                      X_min_fixed,...
                                                      X_max_fixed,...
                                                      options_PSO_2); 
    
                if fval<fval_min_Rx
                    fval_min_Rx = fval;
                    X_out_Rx = X_final(symmetries_all_to_pass);
                    options_PSO_2.SwarmSize = min(options_PSO_2.SwarmSize+100,1000);
                    options_PSO_2.MaxIterations = min(options_PSO_2.MaxIterations+100,1500);
                    options_PSO_2.InertiaRange = [max(0.3,options_PSO_2.InertiaRange(1)-0.1) min(options_PSO_2.InertiaRange(2)+0.1,2.0)];
                    options_PSO_2.SelfAdjustmentWeight = max(0.2,options_PSO_2.SelfAdjustmentWeight-0.1);
                    options_PSO_2.SocialAdjustmentWeight = min(options_PSO_2.SocialAdjustmentWeight+0.1,2.5);
                    Xmin_raw = low_lv * X_final;
                    Xmax_raw = upp_lv * X_final;
                    X_min_fixed = max(min(min(Xmin_raw, Xmax_raw), 1), -1);
                    X_max_fixed = max(min(max(Xmin_raw, Xmax_raw), 1), -1);
                end  
    
                if fval_min_Rx < 1
                    break
                end     
    
            end
    
        else
            X_out_Rx = elems_after_matching(symmetries_all);
            X_out_Rx = rmmissing(X_out_Rx(Rx_elems_detailed));
        end
    end

    if Tx_flag
        Tx_elems_match                                          = match_split==2;
        Tx_elems_tune                                           = tune_split==2;
        Tx_elems_matching_mask_split                            = matching_mask_split==2;
        Tx_elems_fixed_matching_mask_split                      = fixed_matching_mask_split==2;
        Tx_elems_tuning_mask_split                              = tuning_mask_split==2;
        Tx_elems_mutual_mask_split                              = mutual_mask_split==2;
        symmetries_all_to_pass                                  = symmetries_all(Tx_elems_detailed)';
        [~, ~, symmetries_all_to_pass]                          = unique(symmetries_all_to_pass, 'stable'); 

        fval_match = sum(fval_match_all(Tx_elems_match));
        if fval_match > 1
            low_lv = 0.5;
            upp_lv = 1.5;
        else
            low_lv = 0.8;
            upp_lv = 1.2;
        end

        if ~isempty(nonzeros(Tx_elems_detailed))
            
            [y_ids,~]   = find(mask_split == 2);
            non_ids     = setdiff(1:size(YPn, 1), y_ids);
            YPn_to_pass = YPn(y_ids,y_ids)-YPn(y_ids,non_ids)/YPn(non_ids,non_ids)*YPn(non_ids,y_ids);

            M                                  = nonzeros(all_matching_mask(Tx_elems_matching_mask_split));
            F                                  = nonzeros(all_fixed_matching_mask(Tx_elems_fixed_matching_mask_split));
            N                                  = nonzeros(all_tuning_mask(Tx_elems_tuning_mask_split));
            D1                                 = nonzeros(all_mutual_mask(Tx_elems_mutual_mask_split,1));
            D2                                 = nonzeros(all_mutual_mask(Tx_elems_mutual_mask_split,2));
            D1                                 = ismember(N,D1);
            D2                                 = ismember(N,D2);
            [M,F,N]                            = map_to_local(M,F,N);
            FM                                 = sort([M(:); F(:)]);
            [M_new,~]                          = map_to_local_2(M,F);
            D1                                 = N(D1);
            D2                                 = N(D2);
            all_values                         = [D1 D2];                     
            [~, sort_idx]                      = sort(all_values);  
            ranks                              = zeros(size(all_values));           
            ranks(sort_idx)                    = 1:length(all_values); 
            mapped_D1                          = ranks(1:length(D1));
            mapped_D2                          = ranks(length(D1)+1:end);
            length_matching_network            = unique(max(all_port_order(:,Tx_elems_match)));
            idx12                              = sub2ind(size(YPn), D1, D2);
            idx21                              = sub2ind(size(YPn), D2, D1);
            M_detailed                         = nonzeros(all_matching_mask_detailed(Tx_elems_detailed));
            F_detailed                         = nonzeros(all_fixed_matching_mask_detailed(Tx_elems_detailed));
            N_detailed                         = nonzeros(all_tuning_mask_detailed(Tx_elems_detailed));
            [M_detailed,F_detailed,N_detailed] = map_to_local(M_detailed,F_detailed,N_detailed);
            MN_detailed                        = sort([M_detailed(:); N_detailed(:)]);
            MF_detailed                        = sort([M_detailed(:); F_detailed(:)]);
            symmetries_all_to_pass(F_detailed) = [];
            [unique_vals, ~, new_indices]      = unique(nonzeros(symmetries_all_to_pass), 'stable'); 
            new_mapping                        = 1:length(unique_vals);
            symmetries_all_to_pass             = new_mapping(new_indices);

            costFun_match_tune_decoupling = @(X) matching_and_tuning_and_decoupling(X,...
                                                                                    YPn_to_pass,...
                                                                                    symmetries_all_to_pass,...
                                                                                    M_detailed,...
                                                                                    F_detailed,...
                                                                                    N_detailed,...
                                                                                    MN_detailed,...
                                                                                    MF_detailed,...
                                                                                    all_fixed_matching_values(Tx_elems_detailed),...
                                                                                    FM,...
                                                                                    M_new,...
                                                                                    N,...
                                                                                    idx12,...
                                                                                    idx21,...
                                                                                    mapped_D1,...
                                                                                    mapped_D2,...
                                                                                    length_matching_network,...
                                                                                    all_matching_load_string(:,Tx_elems_match),...
                                                                                    all_tuning_load_string(Tx_elems_tune),...
                                                                                    all_port_order(:,Tx_elems_match),...
                                                                                    emc.omega,...
                                                                                    emc.Z0);
            
            fval_min_Tx = inf;
            if fval_match > 50
                lb_t = lower_bound_mat(unique_indices_all);
                ub_t = upper_bound_mat(unique_indices_all);
                lub_m = Tx_elems_detailed(unique_indices_all);
                X_min_fixed = rmmissing(lb_t(lub_m));
                X_max_fixed = rmmissing(ub_t(lub_m));
            else
                lub_m = Tx_elems_detailed(unique_indices_all);
                Xmin_raw = low_lv * rmmissing(elems_after_matching(lub_m));
                Xmax_raw = upp_lv * rmmissing(elems_after_matching(lub_m));
                X_min_fixed = max(min(min(Xmin_raw, Xmax_raw), 1), -1);
                X_max_fixed = max(min(max(Xmin_raw, Xmax_raw), 1), -1);
            end
    
            for counter_repeat = 1:10
                if isempty(X_min_fixed)
                    X_out_Tx = [];
                    break;
                end
                if fval_min_Tx < 1 
                    break
                end
                rng(counter_repeat);
    
                [X_final, fval, ~, ~] = particleswarm(costFun_match_tune_decoupling,...
                                                      length(X_min_fixed),...
                                                      X_min_fixed,...
                                                      X_max_fixed,...
                                                      options_PSO_2);    
                
                if fval<fval_min_Tx
                    fval_min_Tx = fval;
                    X_out_Tx = X_final(symmetries_all_to_pass);
                    options_PSO_2.SwarmSize = min(options_PSO_2.SwarmSize+100,1000);
                    options_PSO_2.MaxIterations = min(options_PSO_2.MaxIterations+100,1500);
                    options_PSO_2.InertiaRange = [max(0.3,options_PSO_2.InertiaRange(1)-0.1) min(options_PSO_2.InertiaRange(2)+0.1,2.0)];
                    options_PSO_2.SelfAdjustmentWeight = max(0.2,options_PSO_2.SelfAdjustmentWeight-0.1);
                    options_PSO_2.SocialAdjustmentWeight = min(options_PSO_2.SocialAdjustmentWeight+0.1,2.5);
                    Xmin_raw = low_lv * X_final;
                    Xmax_raw = upp_lv * X_final;
                    X_min_fixed = max(min(min(Xmin_raw, Xmax_raw), 1), -1);
                    X_max_fixed = max(min(max(Xmin_raw, Xmax_raw), 1), -1);
                end  
                
                if fval_min_Tx < 1
                    break
                end  
    
            end
    
        else
            X_out_Tx = elems_after_matching(symmetries_all);
            X_out_Tx = rmmissing(X_out_Tx(Tx_elems_detailed));
        end
    end

    if TxRx_flag
        TxRx_elems_match                                          = match_split==3;
        TxRx_elems_tune                                           = tune_split==3;
        TxRx_elems_matching_mask_split                            = matching_mask_split==3;
        TxRx_elems_fixed_matching_mask_split                      = fixed_matching_mask_split==3;
        TxRx_elems_tuning_mask_split                              = tuning_mask_split==3;
        TxRx_elems_mutual_mask_split                              = mutual_mask_split==3;
        symmetries_all_to_pass                                    = symmetries_all(TxRx_elems_detailed)';
        [~, ~, symmetries_all_to_pass]                            = unique(symmetries_all_to_pass, 'stable'); 

        fval_match = sum(fval_match_all(TxRx_elems_match));
        if fval_match > 1
            low_lv = 0.5;
            upp_lv = 1.5;
        else
            low_lv = 0.8;
            upp_lv = 1.2;
        end

        if sum(TxRx_elems_match)>1
            
            [y_ids,~]   = find(mask_split == 3);
            non_ids     = setdiff(1:size(YPn, 1), y_ids);
            YPn_to_pass = YPn(y_ids,y_ids)-YPn(y_ids,non_ids)/YPn(non_ids,non_ids)*YPn(non_ids,y_ids);

            M                                  = nonzeros(all_matching_mask(TxRx_elems_matching_mask_split));
            F                                  = nonzeros(all_fixed_matching_mask(TxRx_elems_fixed_matching_mask_split));
            N                                  = nonzeros(all_tuning_mask(TxRx_elems_tuning_mask_split));
            D1                                 = nonzeros(all_mutual_mask(TxRx_elems_mutual_mask_split,1));
            D2                                 = nonzeros(all_mutual_mask(TxRx_elems_mutual_mask_split,2));
            D1                                 = ismember(N,D1);
            D2                                 = ismember(N,D2);
            [M,F,N]                            = map_to_local(M,F,N);
            FM                                 = sort([M(:); F(:)]);
            [M_new,~]                          = map_to_local_2(M,F);
            D1                                 = N(D1);
            D2                                 = N(D2);
            all_values                         = [D1 D2];                     
            [~, sort_idx]                      = sort(all_values);  
            ranks                              = zeros(size(all_values));           
            ranks(sort_idx)                    = 1:length(all_values); 
            mapped_D1                          = ranks(1:length(D1));
            mapped_D2                          = ranks(length(D1)+1:end);
            length_matching_network            = unique(max(all_port_order(:,TxRx_elems_match)));
            idx12                              = sub2ind(size(YPn), D1, D2);
            idx21                              = sub2ind(size(YPn), D2, D1);
            M_detailed                         = nonzeros(all_matching_mask_detailed(TxRx_elems_detailed));
            F_detailed                         = nonzeros(all_fixed_matching_mask_detailed(TxRx_elems_detailed));
            N_detailed                         = nonzeros(all_tuning_mask_detailed(TxRx_elems_detailed));
            [M_detailed,F_detailed,N_detailed] = map_to_local(M_detailed,F_detailed,N_detailed);
            MN_detailed                        = sort([M_detailed(:); N_detailed(:)]);
            MF_detailed                        = sort([M_detailed(:); F_detailed(:)]);
            symmetries_all_to_pass(F_detailed) = [];
            [unique_vals, ~, new_indices]      = unique(nonzeros(symmetries_all_to_pass), 'stable'); 
            new_mapping                        = 1:length(unique_vals);
            symmetries_all_to_pass             = new_mapping(new_indices);

            costFun_match_tune_preamplifierdecoupling = @(X) matching_and_tuning_and_decoupling_and_preampdecoupling(X,...
                                                                                                                     YPn_to_pass,...
                                                                                                                     symmetries_all_to_pass,...
                                                                                                                     M_detailed,...
                                                                                                                     F_detailed,...
                                                                                                                     N_detailed,...
                                                                                                                     MN_detailed,...
                                                                                                                     MF_detailed,...
                                                                                                                     all_fixed_matching_values(TxRx_elems_detailed),...
                                                                                                                     FM,...
                                                                                                                     M_new,...
                                                                                                                     N,...
                                                                                                                     idx12,...
                                                                                                                     idx21,...
                                                                                                                     mapped_D1,...
                                                                                                                     mapped_D2,...
                                                                                                                     length_matching_network,...
                                                                                                                     all_matching_load_string(:,TxRx_elems_match),...
                                                                                                                     all_tuning_load_string(TxRx_elems_tune),...
                                                                                                                     all_port_order(:,TxRx_elems_match),...
                                                                                                                     emc.omega,...
                                                                                                                     emc.Z0,...
                                                                                                                     emc.Preamp_res);
 
            fval_min_TxRx = inf;
            if fval_match > 50
                lb_t = lower_bound_mat(unique_indices_all);
                ub_t = upper_bound_mat(unique_indices_all);
                lub_m = TxRx_elems_detailed(unique_indices_all);
                X_min_fixed = rmmissing(lb_t(lub_m));
                X_max_fixed = rmmissing(ub_t(lub_m));
            else
                lub_m = TxRx_elems_detailed(unique_indices_all);
                Xmin_raw = low_lv * rmmissing(elems_after_matching(lub_m));
                Xmax_raw = upp_lv * rmmissing(elems_after_matching(lub_m));
                X_min_fixed = max(min(min(Xmin_raw, Xmax_raw), 1), -1);
                X_max_fixed = max(min(max(Xmin_raw, Xmax_raw), 1), -1);
            end
    
            for counter_repeat = 1:10
                if isempty(X_min_fixed)
                    X_out_TxRx = [];
                    break;
                end
                if fval_min_TxRx < 1
                    break
                end
                rng(counter_repeat);
    
                [X_final, fval, ~, ~] = particleswarm(costFun_match_tune_preamplifierdecoupling,...
                                                      length(X_min_fixed),...
                                                      X_min_fixed,...
                                                      X_max_fixed,...
                                                      options_PSO_2);  
    
                if fval<fval_min_TxRx
                    fval_min_TxRx = fval;
                    X_out_TxRx = X_final(symmetries_all_to_pass);
                    options_PSO_2.SwarmSize = min(options_PSO_2.SwarmSize+100,1000);
                    options_PSO_2.MaxIterations = min(options_PSO_2.MaxIterations+100,1500);
                    options_PSO_2.InertiaRange = [max(0.3,options_PSO_2.InertiaRange(1)-0.1) min(options_PSO_2.InertiaRange(2)+0.1,2.0)];
                    options_PSO_2.SelfAdjustmentWeight = max(0.2,options_PSO_2.SelfAdjustmentWeight-0.1);
                    options_PSO_2.SocialAdjustmentWeight = min(options_PSO_2.SocialAdjustmentWeight+0.1,2.5);
                    Xmin_raw = low_lv * X_final;
                    Xmax_raw = upp_lv * X_final;
                    X_min_fixed = max(min(min(Xmin_raw, Xmax_raw), 1), -1);
                    X_max_fixed = max(min(max(Xmin_raw, Xmax_raw), 1), -1);
                end  
    
                if fval_min_TxRx < 1
                    break
                end               
            end
    
        else
            X_out_TxRx = elems_after_matching(symmetries_all);
            X_out_TxRx = rmmissing(X_out_TxRx(TxRx_elems_detailed));
        end
    end

    X_out = zeros(length(symmetries_all),1);
    fval_min = 0;
    c = 0;
    if Rx_flag && ~isempty(X_out_Rx)
        X_out(Rx_elems_detailed.*fixed_mask_detailed_split>0) = X_out_Rx;
        c = c + 1;
        if sum(Rx_elems_match)>1
            fval_min = fval_min+fval_min_Rx;
        else
            fval_min = fval_min+sum(fval_match_all(Rx_elems_match));
        end
    end
    if Tx_flag && ~isempty(X_out_Tx)
        X_out(Tx_elems_detailed.*fixed_mask_detailed_split>0) = X_out_Tx;
        c = c + 1;
        if sum(Tx_elems_match)>1
            fval_min = fval_min+fval_min_Tx;
        else
            fval_min = fval_min+sum(fval_match_all(Tx_elems_match));
        end
    end
    if TxRx_flag && ~isempty(X_out_TxRx)
        X_out(TxRx_elems_detailed.*fixed_mask_detailed_split>0) = X_out_TxRx;
        c = c + 1;
        if sum(TxRx_elems_match)>1
            fval_min = fval_min+fval_min_TxRx;
        else
            fval_min = fval_min+sum(fval_match_all(TxRx_elems_match));
        end
    end
    X_out = nonzeros(X_out(unique_indices_all));
    fval_min = fval_min/c;

    if exist('X_out_Rx')
        check_Rx = ~isempty(X_out_Rx);
    else
        check_Rx = 0;
    end
    if exist('X_out_Tx')
        check_Tx = ~isempty(X_out_Tx);
    else
        check_Tx = 0;
    end
    if exist('X_out_TxRx')
        check_TxRx = ~isempty(X_out_TxRx);
    else
        check_TxRx = 0;
    end

end
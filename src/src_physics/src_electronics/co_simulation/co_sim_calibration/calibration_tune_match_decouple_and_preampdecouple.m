function [M_cal,...
          YPm,...
          ZPm,...
          SPm,...
          SP_check,...
          phi_lumped_elements] = calibration_tune_match_decouple_and_preampdecouple(X,...
                                                                                    YPn,...
                                                                                    symmetries,...
                                                                                    matching_mask_detailed,...
                                                                                    fixed_matching_mask_detailed,...
                                                                                    fixed_matching_elements,...
                                                                                    tuning_mask_detailed,...
                                                                                    matching_mask,...
                                                                                    fixed_matching_mask,...
                                                                                    tuning_mask,...
                                                                                    mutual_mask,...
                                                                                    matching_load_string,...
                                                                                    tuning_load_string,...
                                                                                    port_order,...
                                                                                    omega,...
                                                                                    z0,...
                                                                                    preamp_res,...
                                                                                    Q_val)
      
    M                                  = nonzeros(matching_mask);
    F                                  = nonzeros(fixed_matching_mask);
    N                                  = nonzeros(tuning_mask);
    D1                                 = nonzeros(mutual_mask(:,1));
    D2                                 = nonzeros(mutual_mask(:,2));
    D1                                 = ismember(N,D1);
    D2                                 = ismember(N,D2);
    [M,F,N]                            = map_to_local(M,F,N);
    M                                  = sort([M(:); F(:)]);
    D1                                 = N(D1);
    D2                                 = N(D2);
    all_values                         = [D1 D2];                     
    [~, sort_idx]                      = sort(all_values);  
    ranks                              = zeros(size(all_values));           
    ranks(sort_idx)                    = 1:length(all_values); 
    mapped_D1                          = ranks(1:length(D1));
    mapped_D2                          = ranks(length(D1)+1:end);
    length_matching_network            = unique(max(port_order));
    M_detailed                         = nonzeros(matching_mask_detailed);
    F_detailed                         = nonzeros(fixed_matching_mask_detailed);
    N_detailed                         = nonzeros(tuning_mask_detailed);
    [M_detailed,F_detailed,N_detailed] = map_to_local(M_detailed,F_detailed,N_detailed);
    MN_detailed                        = sort([M_detailed(:); N_detailed(:)]);
    MF_detailed                        = sort([M_detailed(:); F_detailed(:)]);

    lumped_elements                  = X(symmetries); 
    all_elements                     = zeros(numel(M_detailed)+numel(F_detailed)+numel(N_detailed),1);
    all_elements(MN_detailed)        = lumped_elements;
    all_elements(F_detailed)         = fixed_matching_elements(F_detailed);
    matching_elements                = all_elements(MF_detailed);
    Q_matching_elements              = Q_val(MF_detailed);
    true_tuning_elements             = all_elements(N_detailed);
    Q_true_tuning_elements           = Q_val(N_detailed);
    is_capacitor_parallel            = strcmpi(matching_load_string, 'capacitorParallel');
    is_inductor_parallel             = strcmpi(matching_load_string, 'inductorParallel');
    is_resistor_parallel             = strcmpi(matching_load_string, 'resistorParallel');
    is_capacitor_series              = strcmpi(matching_load_string, 'capacitorSeries');
    is_inductor_series               = strcmpi(matching_load_string, 'inductorSeries');
    is_resistor_series               = strcmpi(matching_load_string, 'resistorSeries');
    is_capacitor                     = strcmpi(tuning_load_string, 'capacitor');
    is_inductor                      = strcmpi(tuning_load_string, 'inductor');
    is_resistor                      = strcmpi(tuning_load_string, 'resistor');
    is_mutual_inductor               = strcmpi(tuning_load_string, 'mutual_inductor');
    is_mutual_inductance_coefficient = strcmpi(tuning_load_string, 'mutual_inductance_coefficient');
    mutual_inductance_coefficients   = lumped_elements(end-sum(is_mutual_inductance_coefficient)+1:end);
    mutual_inductances               = true_tuning_elements(is_mutual_inductor);

    % Tuning Elements
    Y_vals               = zeros(size(true_tuning_elements));
    Y_vals_loss          = zeros(size(true_tuning_elements));

    Cap_loss             = omega .* true_tuning_elements(is_capacitor) .* Q_true_tuning_elements(is_capacitor);
    Ind_loss             = (omega .* true_tuning_elements(is_inductor) .* Q_true_tuning_elements(is_inductor) );
    MutInd_loss          = (omega .* true_tuning_elements(is_mutual_inductor) .* Q_true_tuning_elements(is_mutual_inductor) );

    Cap_Value            = 1i * omega .* true_tuning_elements(is_capacitor);
    Ind_Value            = 1 ./ (1i * omega .* true_tuning_elements(is_inductor));
    MutInd_Value         = 1 ./ (1i*omega .* true_tuning_elements(is_mutual_inductor));
    
    Y_vals(is_capacitor) = Cap_Value .* Cap_loss ./ (Cap_Value + Cap_loss);
    Y_vals(is_inductor)  = Ind_Value .* Ind_loss ./ (Ind_Value + Ind_loss);
    Y_vals(is_resistor)  = 1 ./ true_tuning_elements(is_resistor);
    Y_vals(is_mutual_inductor) = MutInd_Value .* MutInd_loss ./ (MutInd_Value + MutInd_loss);

    Y_vals_loss(is_capacitor) = Cap_Value;
    Y_vals_loss(is_inductor)  = Ind_Value;
    Y_vals_loss(is_mutual_inductor) = MutInd_Value;

    YPn_loss = YPn;
    
    % Construct diagonal admittance matrix and add values
    E_tu     = diag(Y_vals(:));
    YPn(N,N) = YPn(N,N) + E_tu;
    idx      = sub2ind(size(YPn), D1, D2);
    YPn(idx) = YPn(idx) + 1 ./ (1i*omega .* mutual_inductance_coefficients.*sqrt(mutual_inductances(mapped_D1).*mutual_inductances(mapped_D2)));
    idx      = sub2ind(size(YPn), D2, D1);
    YPn(idx) = YPn(idx) + 1 ./ (1i*omega .* mutual_inductance_coefficients.*sqrt(mutual_inductances(mapped_D2).*mutual_inductances(mapped_D1)));
    ZPn      = inv(YPn);

    M_tune                = zeros(length(N)+length(M),length(M));
    M_tune(M,1:length(M)) = ZPn(M,M);
    M_tune(N,1:length(M)) = ZPn(N,M);
    M_tune                = M_tune/ZPn(M,M);

    % Add tuning elements
    YPm = YPn(M,M)-YPn(M,N)/YPn(N,N)*YPn(N,M);

    % Construct diagonal loss admittance matrix and add values
    E_tu_loss     = diag(Y_vals_loss(:));
    YPn_loss(N,N) = YPn_loss(N,N) + E_tu_loss;
    idx      = sub2ind(size(YPn_loss), D1, D2);
    YPn_loss(idx) = YPn_loss(idx) + 1 ./ (1i*omega .* mutual_inductance_coefficients.*sqrt(mutual_inductances(mapped_D1).*mutual_inductances(mapped_D2)));
    idx      = sub2ind(size(YPn_loss), D2, D1);
    YPn_loss(idx) = YPn_loss(idx) + 1 ./ (1i*omega .* mutual_inductance_coefficients.*sqrt(mutual_inductances(mapped_D2).*mutual_inductances(mapped_D1)));

    % Add tuning elements to loss
    YPm_loss = YPn_loss(M,M)-YPn_loss(M,N)/YPn_loss(N,N)*YPn_loss(N,M);

    % Model preamplifier effect
    P             = size(YPm,1);
    SP_check      = zeros(P,1);
    ZP_check      = zeros(P,1);
    ZP_check_loss = zeros(P,1);
    M_Rx_match    = zeros(P,P);
    M_Rx_decouple = zeros(P,P);
    for k = 1:P
        rps = setdiff(1:P, k);
        YPn_temp = YPm;
        YPn_tempQ = YPm_loss;
         % Decouple
        YPn_temp(rps, rps) = YPn_temp(rps, rps) + 1/preamp_res;
        YPn_tempQ(rps, rps) = YPn_tempQ(rps, rps) + 1/preamp_res;
        ZPn_temp = inv(YPn_temp);
        decouple_col = zeros(P,1);
        decouple_col(k,1) = ZPn_temp(k,k);
        decouple_col(rps,1) = ZPn_temp(rps,k);
        decouple_col = decouple_col / ZPn_temp(k,k);
        M_Rx_decouple(:, k) = decouple_col;
        % Matching Elements
        YPm_temp = YPn_temp(k,k)-YPn_temp(k,rps)/YPn_temp(rps,rps)*YPn_temp(rps,k);
        YPm_tempQ = YPn_tempQ(k,k)-YPn_tempQ(k,rps)/YPn_tempQ(rps,rps)*YPn_tempQ(rps,k);
        SPs = np_z2s(inv(YPm_temp),z0);
        for i = 1:length_matching_network
            stage_mask            = (port_order == i);
            existsXx              = is_capacitor_parallel | is_inductor_parallel | is_resistor_parallel | is_capacitor_series | is_inductor_series | is_resistor_series;
            slot2idx              = zeros(numel(existsXx),1);
            slot2idx(existsXx)    = 1:nnz(existsXx);          
            ii                    = slot2idx(stage_mask); 
            stage_elements        = zeros(nnz(stage_mask),1);
            stage_elementsQ       = zeros(nnz(stage_mask),1);
            good                  = ii > 0;
            stage_elements(good)  = matching_elements(ii(good));
            stage_elementsQ(good) = Q_matching_elements(ii(good));
            is_capacitor_parallel_local = is_capacitor_parallel(port_order==i);
            is_inductor_parallel_local  = is_inductor_parallel(port_order==i);
            is_resistor_parallel_local  = is_resistor_parallel(port_order==i);
            is_capacitor_series_local   = is_capacitor_series(port_order==i);
            is_inductor_series_local    = is_inductor_series(port_order==i);
            is_resistor_series_local    = is_resistor_series(port_order==i);
            if sum(is_capacitor_parallel_local)> 0 || sum(is_inductor_parallel_local)>0 || sum(is_resistor_parallel_local)>0
                E_cp = diag(1i* omega * stage_elements(is_capacitor_parallel_local)); 
                E_lp = diag(1./(1i*omega*stage_elements(is_inductor_parallel_local)));    
                E_rp = diag(1./(stage_elements(is_resistor_parallel_local)));    
            
                Cap_loss  = diag(omega * stage_elements(is_capacitor_parallel_local) .* stage_elementsQ(is_capacitor_parallel_local) );
                Ind_loss  = diag((omega * stage_elements(is_inductor_parallel_local).* stage_elementsQ(is_inductor_parallel_local) ));    

                if ~isempty(E_cp)
                    YPm_temp = YPm_temp + E_cp(k,k).*Cap_loss(k,k)./(E_cp(k,k)+Cap_loss(k,k));
                    YPm_tempQ = YPm_tempQ + E_cp(k,k);
                end
                if ~isempty(E_lp)
                    YPm_temp = YPm_temp + E_lp(k,k).*Ind_loss(k,k)./(E_lp(k,k)+Ind_loss(k,k));
                    YPm_tempQ = YPm_tempQ + E_lp(k,k);
                end
                if ~isempty(E_rp)
                    YPm_temp = YPm_temp + E_rp(k,k);
                end
            end
            if sum(is_capacitor_series_local)> 0 || sum(is_inductor_series_local)>0 || sum(is_resistor_series_local)>0
                E_cs = diag(1./(1i*omega*stage_elements(is_capacitor_series_local))); 
                E_ls = diag(1i*omega*stage_elements(is_inductor_series_local));    
                E_rs = diag(stage_elements(is_resistor_series_local));    
            
                Cap_loss  = diag((omega * stage_elements(is_capacitor_series_local).* stage_elementsQ(is_capacitor_series_local) ));   
                Ind_loss  = diag(omega * stage_elements(is_inductor_series_local) .* stage_elementsQ(is_inductor_series_local) ); 

                ZPm_temp = inv(YPm_temp);
                ZPm_tempQ = inv(YPm_tempQ);
                if ~isempty(E_cs)
                    ZPm_temp = ZPm_temp + E_cs(k,k).*Cap_loss(k,k)./(E_cs(k,k)+Cap_loss(k,k));
                    ZPm_tempQ = ZPm_tempQ + E_cs(k,k);
                end
                if ~isempty(E_ls)
                    ZPm_temp = ZPm_temp + E_ls(k,k).*Ind_loss(k,k)./(E_ls(k,k)+Ind_loss(k,k));
                    ZPm_tempQ = ZPm_tempQ + E_ls(k,k);
                end
                if ~isempty(E_rs)
                    ZPm_temp = ZPm_temp + E_rs(k,k);
                end
                YPm_temp = inv(ZPm_temp);
                YPm_tempQ = inv(ZPm_tempQ);
            end
        end
        SPm = np_z2s(inv(YPm_temp),z0);
        M_Rx_match(k,k) = calibration_matching(SPs,SPm,z0);
        SP_check(k,1) = SPm;
        ZP_check_loss(k,1) = inv(YPm_tempQ);
        ZP_check(k,1) = inv(YPm_temp);
    end

    % Model matching elements
    SPs = np_z2s(inv(YPm),z0);
    for i = 1:length_matching_network
        stage_mask            = (port_order == i);
        existsXx              = is_capacitor_parallel | is_inductor_parallel | is_resistor_parallel | is_capacitor_series | is_inductor_series | is_resistor_series;
        slot2idx              = zeros(numel(existsXx),1);
        slot2idx(existsXx)    = 1:nnz(existsXx);          
        ii                    = slot2idx(stage_mask); 
        stage_elements        = zeros(nnz(stage_mask),1);
        stage_elementsQ       = zeros(nnz(stage_mask),1);
        good                  = ii > 0;
        stage_elements(good)  = matching_elements(ii(good));
        stage_elementsQ(good) = Q_matching_elements(ii(good));
        is_capacitor_parallel_local = is_capacitor_parallel(port_order==i);
        is_inductor_parallel_local  = is_inductor_parallel(port_order==i);
        is_resistor_parallel_local  = is_resistor_parallel(port_order==i);
        is_capacitor_series_local   = is_capacitor_series(port_order==i);
        is_inductor_series_local    = is_inductor_series(port_order==i);
        is_resistor_series_local    = is_resistor_series(port_order==i);
        if sum(is_capacitor_parallel_local)> 0 || sum(is_inductor_parallel_local)>0 || sum(is_resistor_parallel_local)>0
            E_cp = (1i* omega * stage_elements(is_capacitor_parallel_local)); 
            E_lp = (1./(1i*omega*stage_elements(is_inductor_parallel_local)));    
            E_rp = diag(1./(stage_elements(is_resistor_parallel_local)));  
            
            Cap_loss  = (omega * stage_elements(is_capacitor_parallel_local) .* stage_elementsQ(is_capacitor_parallel_local) );
            Ind_loss  = ((omega * stage_elements(is_inductor_parallel_local).* stage_elementsQ(is_inductor_parallel_local) ));  

            if ~isempty(E_cp)
                YPm = YPm + diag(E_cp.*Cap_loss./(E_cp+Cap_loss));
                YPm_loss = YPm_loss + diag(E_cp);
            end
            if ~isempty(E_lp)
                YPm = YPm + diag(E_lp.*Ind_loss./(E_lp+Ind_loss));
                YPm_loss = YPm_loss + diag(E_lp);
            end
            if ~isempty(E_rp)
                YPm = YPm + E_rp;
            end
        end
        if sum(is_capacitor_series_local)> 0 || sum(is_inductor_series_local)>0 || sum(is_resistor_series_local)>0
            E_cs = (1./(1i*omega*stage_elements(is_capacitor_series_local))); 
            E_ls = (1i*omega*stage_elements(is_inductor_series_local));    
            E_rs = diag(stage_elements(is_resistor_series_local));    
            
            Cap_loss  = ((omega * stage_elements(is_capacitor_series_local).* stage_elementsQ(is_capacitor_series_local) ));   
            Ind_loss  = (omega * stage_elements(is_inductor_series_local) .* stage_elementsQ(is_inductor_series_local) ); 

            ZPm = inv(YPm);
            ZPm_loss = inv(YPm_loss);
            if ~isempty(E_cs)
                ZPm = ZPm + diag(E_cs.*Cap_loss./(E_cs+Cap_loss));
                ZPm_loss = ZPm_loss + diag(E_cs);
            end
            if ~isempty(E_ls)
                ZPm = ZPm + diag(E_ls.*Ind_loss./(E_ls+Ind_loss));
                ZPm_loss = ZPm_loss + diag(E_ls);
            end
            if ~isempty(E_rs)
                ZPm = ZPm + E_rs;
            end
            YPm = inv(ZPm);
            YPm_loss = inv(ZPm_loss);
        end
    end
    ZPm = inv(YPm);
    SPm = np_z2s(ZPm,z0);
    M_Tx_match = calibration_matching(SPs,SPm,z0);

    M_cal.rx = (M_tune*M_Rx_decouple)*M_Rx_match;
    M_cal.tx = M_tune*M_Tx_match;

    ZPm_loss = inv(YPm_loss);
    phi_lumped_elements.Tx = diag(1/2*diag(abs(real(ZPm_loss-ZPm))).*(1./abs(diag(ZPm))).^2);
    phi_lumped_elements.Rx = 1/2*diag(abs(real(ZP_check_loss-ZP_check)).*(1./abs(ZP_check)).^2);

end
function [M_cal,...
          YPm_Tx,...
          ZPm,...
          SPm,...
          SP_check,...
          phi_lumped_elements] = calibration_tune_match_Tx_decouple_and_TxRx_preampdecouple(X,...
                                                                                            YPn,...
                                                                                            symmetries,...
                                                                                            match_split,...
                                                                                            mask_detailed_split,...
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
                                                                                            detune_res,...
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
    length_matching_network            = unique(max(port_order(:)));
    M_detailed                         = nonzeros(matching_mask_detailed);
    F_detailed                         = nonzeros(fixed_matching_mask_detailed);
    N_detailed                         = nonzeros(tuning_mask_detailed);
    [M_detailed,F_detailed,N_detailed] = map_to_local(M_detailed,F_detailed,N_detailed);
    MN_detailed                        = sort([M_detailed(:); N_detailed(:)]);
    MF_detailed                        = sort([M_detailed(:); F_detailed(:)]);

    lumped_elements                  = X(symmetries); 
    Rx_elems_match                   = mask_detailed_split==3;
    Tx_elems_match                   = mask_detailed_split==2 | mask_detailed_split==3;
    RxM                              = match_split==3;
    TxM_det                          = match_split==2;
    TxM                              = match_split==2 | match_split==3;
    all_elements                     = zeros(numel(M_detailed)+numel(F_detailed)+numel(N_detailed),1);
    all_elements(MN_detailed)        = lumped_elements;
    all_elements(F_detailed)         = nonzeros(fixed_matching_elements);
    true_tuning_elements             = all_elements(N_detailed);
    Q_true_tuning_elements           = Q_val(N_detailed);
    matching_elementsRx              = all_elements(MF_detailed(Rx_elems_match(MF_detailed)));
    matching_elementsTx              = all_elements(MF_detailed(Tx_elems_match(MF_detailed)));
    Q_matching_elementsRx            = Q_val(MF_detailed(Rx_elems_match(MF_detailed)));
    Q_matching_elementsTx            = Q_val(MF_detailed(Tx_elems_match(MF_detailed)));
    matching_load_stringRx           = matching_load_string(:,RxM);
    matching_load_stringTx           = matching_load_string(:,TxM);
    port_orderRx                     = port_order(:,RxM);
    port_orderRx                     = port_orderRx(:);
    port_orderTx                     = port_order(:,TxM);
    port_orderTx                     = port_orderTx(:);
    is_capacitor_parallelRx          = strcmpi(matching_load_stringRx(:), 'capacitorParallel');
    is_inductor_parallelRx           = strcmpi(matching_load_stringRx(:), 'inductorParallel');
    is_resistor_parallelRx           = strcmpi(matching_load_stringRx(:), 'resistorParallel');
    is_capacitor_seriesRx            = strcmpi(matching_load_stringRx(:), 'capacitorSeries');
    is_inductor_seriesRx             = strcmpi(matching_load_stringRx(:), 'inductorSeries');
    is_resistor_seriesRx             = strcmpi(matching_load_stringRx(:), 'resistorSeries');
    is_capacitor_parallelTx          = strcmpi(matching_load_stringTx(:), 'capacitorParallel');
    is_inductor_parallelTx           = strcmpi(matching_load_stringTx(:), 'inductorParallel');
    is_resistor_parallelTx           = strcmpi(matching_load_stringTx(:), 'resistorParallel');
    is_capacitor_seriesTx            = strcmpi(matching_load_stringTx(:), 'capacitorSeries');
    is_inductor_seriesTx             = strcmpi(matching_load_stringTx(:), 'inductorSeries');
    is_resistor_seriesTx             = strcmpi(matching_load_stringTx(:), 'resistorSeries');
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

    % Construct diagonal loss admittance matrix and add values
    E_tu_loss     = diag(Y_vals_loss(:));
    YPn_loss(N,N) = YPn_loss(N,N) + E_tu_loss;
    idx           = sub2ind(size(YPn_loss), D1, D2);
    YPn_loss(idx) = YPn_loss(idx) + 1 ./ (1i*omega .* mutual_inductance_coefficients.*sqrt(mutual_inductances(mapped_D1).*mutual_inductances(mapped_D2)));
    idx           = sub2ind(size(YPn_loss), D2, D1);
    YPn_loss(idx) = YPn_loss(idx) + 1 ./ (1i*omega .* mutual_inductance_coefficients.*sqrt(mutual_inductances(mapped_D2).*mutual_inductances(mapped_D1)));

    M_tune                = zeros(length(N)+length(M),length(M));
    M_tune(M,1:length(M)) = ZPn(M,M);
    M_tune(N,1:length(M)) = ZPn(N,M);
    M_tune                = M_tune/ZPn(M,M);

    % Add tuning elements
    YPm = YPn(M,M)-YPn(M,N)/YPn(N,N)*YPn(N,M);

    % Apply detuning
    Txdet = 1/detune_res*eye(size(nonzeros(TxM_det),1));
    YPm_Rx = YPm;
    YPm_Rx(TxM_det,TxM_det) = YPm_Rx(TxM_det,TxM_det)+Txdet;
    ZPn_Rx = inv(YPm_Rx);
    
    MRx                        = size(nonzeros(RxM),1);
    M_detune_Rx                = zeros(length(M),MRx);
    M_detune_Rx(RxM,1:MRx)     = ZPn_Rx(RxM,RxM);
    M_detune_Rx(TxM_det,1:MRx) = ZPn_Rx(TxM_det,RxM);
    M_detune_Rx                = M_detune_Rx/ZPn_Rx(RxM,RxM);

    % Switch to de-tuned mode
    YPm_Rx = YPm(RxM,RxM)-YPm(RxM,TxM_det)/(YPm(TxM_det,TxM_det)+Txdet)*YPm(TxM_det,RxM);
    YPm_Tx = YPm;

    % Add tuning elements to loss
    YPm_loss = YPn_loss(M,M)-YPn_loss(M,N)/YPn_loss(N,N)*YPn_loss(N,M);
    % Add tuning elements to loss
    YPm_Rx_loss = YPm_loss(RxM,RxM)-YPm_loss(RxM,TxM_det)/(YPm_loss(TxM_det,TxM_det)+Txdet)*YPm_loss(TxM_det,RxM);
    YPm_Tx_loss = YPm_loss;

    % Model preamplifier effect
    P             = size(YPm_Rx,1);
    SP_check      = zeros(P,1);
    ZP_check      = zeros(P,1);
    ZP_check_loss = zeros(P,1);
    M_Rx_match    = zeros(P,P);
    M_Rx_decouple = zeros(P,P);
    for k = 1:P
        rps = setdiff(1:P, k);
        YPn_temp = YPm_Rx;
        YPn_tempQ = YPm_Rx_loss;
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
            stage_mask            = (port_orderRx == i);
            existsXx              = is_capacitor_parallelRx | is_inductor_parallelRx | is_resistor_parallelRx | is_capacitor_seriesRx | is_inductor_seriesRx | is_resistor_seriesRx;
            slot2idx              = zeros(numel(existsXx),1);
            slot2idx(existsXx)    = 1:nnz(existsXx);          
            ii                    = slot2idx(stage_mask); 
            stage_elements        = zeros(nnz(stage_mask),1);
            stage_elementsQ       = zeros(nnz(stage_mask),1);
            good                  = ii > 0;
            stage_elements(good)  = matching_elementsRx(ii(good));
            stage_elementsQ(good) = Q_matching_elementsRx(ii(good));
            is_capacitor_parallel_local = is_capacitor_parallelRx(port_orderRx==i);
            is_inductor_parallel_local  = is_inductor_parallelRx(port_orderRx==i);
            is_resistor_parallel_local  = is_resistor_parallelRx(port_orderRx==i);
            is_capacitor_series_local   = is_capacitor_seriesRx(port_orderRx==i);
            is_inductor_series_local    = is_inductor_seriesRx(port_orderRx==i);
            is_resistor_series_local    = is_resistor_seriesRx(port_orderRx==i);
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
    SPs = np_z2s(inv(YPm_Tx),z0);
    for i = 1:length_matching_network
        stage_mask            = (port_orderTx == i);
        existsXx              = is_capacitor_parallelTx | is_inductor_parallelTx | is_resistor_parallelTx | is_capacitor_seriesTx | is_inductor_seriesTx | is_resistor_seriesTx;
        slot2idx              = zeros(numel(existsXx),1);
        slot2idx(existsXx)    = 1:nnz(existsXx);          
        ii                    = slot2idx(stage_mask); 
        stage_elements        = zeros(nnz(stage_mask),1);
        stage_elementsQ       = zeros(nnz(stage_mask),1);
        good                  = ii > 0;
        stage_elements(good)  = matching_elementsTx(ii(good));
        stage_elementsQ(good) = Q_matching_elementsTx(ii(good));
        is_capacitor_parallel_local = is_capacitor_parallelTx(port_orderTx==i);
        is_inductor_parallel_local  = is_inductor_parallelTx(port_orderTx==i);
        is_resistor_parallel_local  = is_resistor_parallelTx(port_orderTx==i);
        is_capacitor_series_local   = is_capacitor_seriesTx(port_orderTx==i);
        is_inductor_series_local    = is_inductor_seriesTx(port_orderTx==i);
        is_resistor_series_local    = is_resistor_seriesTx(port_orderTx==i);
        if sum(is_capacitor_parallel_local)> 0 || sum(is_inductor_parallel_local)>0 || sum(is_resistor_parallel_local)>0
            E_cp = (1i* omega * stage_elements(is_capacitor_parallel_local)); 
            E_lp = (1./(1i*omega*stage_elements(is_inductor_parallel_local)));    
            E_rp = diag(1./(stage_elements(is_resistor_parallel_local)));    
            
            Cap_loss  = (omega * stage_elements(is_capacitor_parallel_local) .* stage_elementsQ(is_capacitor_parallel_local) );
            Ind_loss  = ((omega * stage_elements(is_inductor_parallel_local).* stage_elementsQ(is_inductor_parallel_local) ));  

            if ~isempty(E_cp)
                YPm_Tx = YPm_Tx + diag(E_cp.*Cap_loss./(E_cp+Cap_loss));
                YPm_Tx_loss = YPm_Tx_loss + diag(E_cp);
            end
            if ~isempty(E_lp)
                YPm = YPm + diag(E_lp.*Ind_loss./(E_lp+Ind_loss));
                YPm_Tx_loss = YPm_Tx_loss + diag(E_lp);
            end
            if ~isempty(E_rp)
                YPm_Tx = YPm_Tx + E_rp;
            end
        end
        if sum(is_capacitor_series_local)> 0 || sum(is_inductor_series_local)>0 || sum(is_resistor_series_local)>0
            E_cs = (1./(1i*omega*stage_elements(is_capacitor_series_local))); 
            E_ls = (1i*omega*stage_elements(is_inductor_series_local));    
            E_rs = diag(stage_elements(is_resistor_series_local));    
            
            Cap_loss  = ((omega * stage_elements(is_capacitor_series_local).* stage_elementsQ(is_capacitor_series_local) ));   
            Ind_loss  = (omega * stage_elements(is_inductor_series_local) .* stage_elementsQ(is_inductor_series_local) ); 

            ZPm_Tx = inv(YPm_Tx);
            ZPm_Tx_loss = inv(YPm_Tx_loss);
            if ~isempty(E_cs)
                ZPm_Tx = ZPm_Tx + diag(E_cs.*Cap_loss./(E_cs+Cap_loss));
                ZPm_Tx_loss = ZPm_Tx_loss + diag(E_cs);
            end
            if ~isempty(E_ls)
                ZPm_Tx = ZPm_Tx + diag(E_ls.*Ind_loss./(E_ls+Ind_loss));
                ZPm_Tx_loss = ZPm_Tx_loss + diag(E_ls);
            end
            if ~isempty(E_rs)
                ZPm_Tx = ZPm_Tx + E_rs;
            end
            YPm_Tx = inv(ZPm_Tx);
            YPm_Tx_loss = inv(ZPm_Tx_loss);
        end
    end
    ZPm = inv(YPm_Tx);
    SPm = np_z2s(ZPm,z0);
    M_Tx_match = calibration_matching(SPs,SPm,z0);

    M_cal.rx = ((M_tune*M_detune_Rx)*M_Rx_decouple)*M_Rx_match;
    M_cal.tx = M_tune*M_Tx_match;

    ZPm_loss = inv(YPm_Tx_loss);
    phi_lumped_elements.Tx = diag(1/2*diag(abs(real(ZPm_loss-ZPm))).*(1./abs(diag(ZPm))).^2);
    phi_lumped_elements.Rx = 1/2*diag(abs(real(ZP_check_loss-ZP_check)).*(1./abs(ZP_check)).^2);

end
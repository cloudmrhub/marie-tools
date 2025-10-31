function [cf] = matching_and_tuning_and_TxRx_decoupling_and_Rx_preampdecoupling(X,...
                                                                                YPn,...
                                                                                symmetries,...
                                                                                M_detailed,...
                                                                                F_detailed,...
                                                                                MF_detailed,...
                                                                                N_detailed,...
                                                                                MN_detailed,...
                                                                                fixed_matching_elements,...
                                                                                M,...
                                                                                Final_M,...
                                                                                N,...
                                                                                idx12,...
                                                                                idx21,...
                                                                                mapped_D1,...
                                                                                mapped_D2,...
                                                                                length_matching_network,...
                                                                                matching_load_string,...
                                                                                tuning_load_string,...
                                                                                port_order,...
                                                                                match_split,...
                                                                                mask_detailed_split,...
                                                                                omega,...
                                                                                z0,...
                                                                                preamp_res,...
                                                                                detune_res,...
                                                                                check_Rx,...
                                                                                check_Tx,...
                                                                                check_TxRx)

    lumped_elements                  = X(symmetries); 
    Rx_elems_match                   = mask_detailed_split==1 | mask_detailed_split==3;
    Tx_elems_match                   = mask_detailed_split==3;
    RxM                              = match_split==1 | match_split==3;
    RxM_det                          = match_split==1;
    TxM                              = match_split==3;
    Final_M_Rx                       = Final_M(RxM(Final_M) == 1);
    Final_M_Tx                       = Final_M(TxM(Final_M) == 1);
    [Final_M_Rx, Final_M_Tx]         = map_to_local_2(Final_M_Rx, Final_M_Tx);
    all_elements                     = zeros(numel(M_detailed)+numel(F_detailed)+numel(N_detailed),1);
    all_elements(MN_detailed)        = lumped_elements;
    all_elements(F_detailed)         = nonzeros(fixed_matching_elements);
    true_tuning_elements             = all_elements(N_detailed);
    matching_elementsRx              = all_elements(MF_detailed(Rx_elems_match(MF_detailed)));
    matching_elementsTx              = all_elements(MF_detailed(Tx_elems_match(MF_detailed)));
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
    Y_vals(is_capacitor) = 1i * omega .* true_tuning_elements(is_capacitor);
    Y_vals(is_inductor)  = 1 ./ (1i * omega .* true_tuning_elements(is_inductor));
    Y_vals(is_resistor)  = 1 ./ true_tuning_elements(is_resistor);
    Y_vals(is_mutual_inductor) = 1 ./ (1i*omega .* true_tuning_elements(is_mutual_inductor));
    
    % Construct diagonal admittance matrix and add values
    E_tu     = diag(Y_vals(:));
    YPn(N,N) = YPn(N,N) + E_tu;
    if ~isempty(idx12)
        YPn(idx12) = YPn(idx12) + 1 ./ (1i*omega .* mutual_inductance_coefficients.*sqrt(mutual_inductances(mapped_D1).*mutual_inductances(mapped_D2)));
        YPn(idx21) = YPn(idx21) + 1 ./ (1i*omega .* mutual_inductance_coefficients.*sqrt(mutual_inductances(mapped_D2).*mutual_inductances(mapped_D1)));
    end

    % Aply Schur and get cost-function
    YPm = YPn(M,M)-YPn(M,N)/YPn(N,N)*YPn(N,M);

    % Apply detuning
    Rxdet = 1/detune_res*eye(size(nonzeros(RxM_det),1));
    YPm_Rx = YPm;
    YPm_Tx = YPm(TxM,TxM)-YPm(TxM,RxM_det)/(YPm(RxM_det,RxM_det)+Rxdet)*YPm(RxM_det,TxM);

    % Model preamplifier effect
    P = size(YPm_Rx,1);
    ZP_check = zeros(P,1);
    for k = 1:P
        rps = setdiff(1:P,k);
        YPm_temp = YPm_Rx;
        YPm_temp(rps,rps) = YPm_temp(rps,rps) + 1/preamp_res;
        % Matching Elements
        for i = 1:length_matching_network
            stage_mask           = (port_orderRx == i);
            existsXx             = is_capacitor_parallelRx | is_inductor_parallelRx | is_resistor_parallelRx | is_capacitor_seriesRx | is_inductor_seriesRx | is_resistor_seriesRx;
            slot2idx             = zeros(numel(existsXx),1);
            slot2idx(existsXx)   = 1:nnz(existsXx);          
            ii                   = slot2idx(stage_mask); 
            stage_elements       = zeros(nnz(stage_mask),1);
            good                 = ii > 0;
            stage_elements(good) = matching_elementsRx(ii(good));
            % stage_elements = matching_elementsRx(stage_mask);
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
                if ~isempty(E_cp)
                    YPm_temp(k,k) = YPm_temp(k,k) + E_cp(k,k);
                end
                if ~isempty(E_lp)
                    YPm_temp(k,k) = YPm_temp(k,k) + E_lp(k,k);
                end
                if ~isempty(E_rp)
                    YPm_temp(k,k) = YPm_temp(k,k) + E_rp(k,k);
                end
            end
            if sum(is_capacitor_series_local)> 0 || sum(is_inductor_series_local)>0 || sum(is_resistor_series_local)>0
                E_cs = diag(1./(1i*omega*stage_elements(is_capacitor_series_local))); 
                E_ls = diag(1i*omega*stage_elements(is_inductor_series_local));    
                E_rs = diag(stage_elements(is_resistor_series_local));    
                ZPm_temp = inv(YPm_temp);
                if ~isempty(E_cs)
                    ZPm_temp(k,k) = ZPm_temp(k,k) + E_cs(k,k);
                end
                if ~isempty(E_ls)
                    ZPm_temp(k,k) = ZPm_temp(k,k) + E_ls(k,k);
                end
                if ~isempty(E_rs)
                    ZPm_temp(k,k) = ZPm_temp(k,k) + E_rs(k,k);
                end
                YPm_temp = inv(ZPm_temp);
            end
        end
        ZP_check(k,1) = inv(YPm_temp(k,k)-YPm_temp(k,rps)/YPm_temp(rps,rps)*YPm_temp(rps,k));
    end

    % Model Decoupling
    for i = 1:length_matching_network
        stage_mask           = (port_orderTx == i);
        existsXx             = is_capacitor_parallelTx | is_inductor_parallelTx | is_resistor_parallelTx | is_capacitor_seriesTx | is_inductor_seriesTx | is_resistor_seriesTx;
        slot2idx             = zeros(numel(existsXx),1);
        slot2idx(existsXx)   = 1:nnz(existsXx);          
        ii                   = slot2idx(stage_mask); 
        stage_elements       = zeros(nnz(stage_mask),1);
        good                 = ii > 0;
        stage_elements(good) = matching_elementsTx(ii(good));
        % stage_elements = matching_elementsTx(stage_mask);
        
        is_capacitor_parallel_local = is_capacitor_parallelTx(port_orderTx==i);
        is_inductor_parallel_local  = is_inductor_parallelTx(port_orderTx==i);
        is_resistor_parallel_local  = is_resistor_parallelTx(port_orderTx==i);
        is_capacitor_series_local   = is_capacitor_seriesTx(port_orderTx==i);
        is_inductor_series_local    = is_inductor_seriesTx(port_orderTx==i);
        is_resistor_series_local    = is_resistor_seriesTx(port_orderTx==i);
        
        if sum(is_capacitor_parallel_local)> 0 || sum(is_inductor_parallel_local)>0 || sum(is_resistor_parallel_local)>0
            E_cp = diag(1i* omega * stage_elements(is_capacitor_parallel_local)); 
            E_lp = diag(1./(1i*omega*stage_elements(is_inductor_parallel_local)));    
            E_rp = diag(1./(stage_elements(is_resistor_parallel_local)));   
            if ~isempty(E_cp)
                YPm_Tx = YPm_Tx + E_cp;
            end
            if ~isempty(E_lp)
                YPm_Tx = YPm_Tx + E_lp;
            end
            if ~isempty(E_rp)
                YPm_Tx = YPm_Tx + E_rp;
            end
        end

        if sum(is_capacitor_series_local)> 0 || sum(is_inductor_series_local)>0 || sum(is_resistor_series_local)>0
            E_cs = diag(1./(1i*omega*stage_elements(is_capacitor_series_local))); 
            E_ls = diag(1i*omega*stage_elements(is_inductor_series_local));    
            E_rs = diag(stage_elements(is_resistor_series_local));    
            ZPm_Tx = inv(YPm_Tx);
            if ~isempty(E_cs)
                ZPm_Tx = ZPm_Tx + E_cs;
            end
            if ~isempty(E_ls)
                ZPm_Tx = ZPm_Tx + E_ls;
            end
            if ~isempty(E_rs)
                ZPm_Tx = ZPm_Tx + E_rs;
            end
            YPm_Tx = inv(ZPm_Tx);
        end

    end
    ZPm_Tx = inv(YPm_Tx);

    w1 = 0; w2 = 0; w3 = 0; w0 = 0;
    W  = toeplitz([w1 w2 w3 w0*ones(1,size(ZPm_Tx,2)-3)],[w1 w2 w3 w0*ones(1,size(ZPm_Tx,2)-3)]);
    W  = W(1:size(ZPm_Tx,2),1:size(ZPm_Tx,2));
    if size(ZPm_Tx,2) > 2
        W(end,1) = w2;
    end
    if size(ZPm_Tx,2) > 3
        W(end-1,1) = w3;
    end
    W = tril(W,-1);

    if check_Rx || check_TxRx
        cf_Rx = norm(real(ZP_check(Final_M_Rx))-z0,'fro') + norm(imag(ZP_check),'fro');
    else
        cf_Rx = 0;
    end
    if check_Tx || check_TxRx
        cf_Tx = norm(diag(real(ZPm_Tx(Final_M_Tx,Final_M_Tx)))-z0,'fro') + norm(diag(imag(ZPm_Tx)),'fro') + norm(W.*ZPm_Tx,'fro');
    else
        cf_Tx = 0;
    end

    cf = cf_Tx + cf_Rx;

end
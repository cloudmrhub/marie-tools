function [M_cal,...
          phi_lumped_elements,...
          YPm_Tx,...
          ZPm,...
          SPm,...
          SP_check] = calibration_match_TxRx_decouple_and_Rx_preampdecouple(matching_elements,...
                                                                            YPm,...
                                                                            match_split,...
                                                                            mask_detailed_split,...
                                                                            matching_load_string,...
                                                                            port_order,...
                                                                            omega,...
                                                                            z0,...
                                                                            preamp_res,...
                                                                            detune_res,...
                                                                            Q_val)
    
    Rx_elems_match                   = mask_detailed_split==1 | mask_detailed_split==3;
    Tx_elems_match                   = mask_detailed_split==3;
    matching_elementsRx              = matching_elements(Rx_elems_match);
    matching_elementsTx              = matching_elements(Tx_elems_match);
    Q_valRx                          = Q_val(Rx_elems_match);
    Q_valTx                          = Q_val(Tx_elems_match);
    RxM                              = match_split==1;
    TxM                              = match_split==3;
    matching_load_stringRx           = matching_load_string(Rx_elems_match);
    matching_load_stringTx           = matching_load_string(Tx_elems_match);
    port_orderRx                     = port_order(Rx_elems_match);
    port_orderTx                     = port_order(Tx_elems_match);
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

    length_matching_network = unique(max(port_order));

    % Apply detuning
    Rxdet = 1/detune_res*eye(size(nonzeros(RxM),1));
    YPm_Tx = YPm;
    YPm_Tx(RxM,RxM) = YPm_Tx(RxM,RxM)+Rxdet;
    ZPn_Tx = inv(YPm_Tx);
    
    MTx                    = size(nonzeros(TxM),1);
    M_detune_Tx            = zeros(length(M),MTx);
    M_detune_Tx(TxM,1:MTx) = ZPn_Tx(TxM,TxM);
    M_detune_Tx(RxM,1:MTx) = ZPn_Tx(RxM,TxM);
    M_detune_Tx            = M_detune_Tx/ZPn_Tx(TxM,TxM);

    % Switch to de-tuned mode
    YPm_Rx = YPm;
    YPm_Tx = YPm(TxM,TxM)-YPm(TxM,RxM)/(YPm(RxM,RxM)+Rxdet)*YPm(RxM,TxM);

    YPm_loss_Rx = YPm_Rx;
    YPm_loss_Tx = YPm_Tx;
    % Model preamplifier effect
    P = size(YPm_Rx,1);
    SP_check = zeros(P,1);
    ZP_check      = zeros(P,1);
    ZP_check_loss = zeros(P,1);
    M_match_Rx = zeros(P,P);
    M_decouple_Rx = zeros(P,P);
    for k = 1:P
        rps = setdiff(1:P, k);
        YPn_temp = YPm_Rx;
        YPn_tempQ = YPm_loss_Rx;
         % Decouple
        YPn_temp(rps, rps) = YPn_temp(rps, rps) + 1/preamp_res;
        YPn_tempQ(rps, rps) = YPn_tempQ(rps, rps) + 1/preamp_res;
        ZPn_temp = inv(YPn_temp);
        decouple_col = zeros(P,1);
        decouple_col(k,1) = ZPn_temp(k,k);
        decouple_col(rps,1) = ZPn_temp(rps,k);
        decouple_col = decouple_col / ZPn_temp(k,k);
        M_decouple_Rx(:, k) = decouple_col;
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
            stage_elementsQ(good) = Q_valRx(ii(good));
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
                Ind_loss  = diag(1./(omega * stage_elements(is_inductor_parallel_local).* stage_elementsQ(is_inductor_parallel_local) ));     
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
            
                Cap_loss  = diag(1./(omega * stage_elements(is_capacitor_series_local).* stage_elementsQ(is_capacitor_series_local) ));   
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
        M_match_Rx(k,k) = calibration_matching(SPs,SPm,z0);
        SP_check(k,1) = SPm;
        ZP_check_loss(k,1) = inv(YPm_tempQ);
        ZP_check(k,1) = inv(YPm_temp);
    end

    SPs = np_z2s(inv(YPm_Tx),z0);
    % Model matching elements
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
        stage_elementsQ(good) = Q_valTx(ii(good));
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
            Ind_loss  = (1./(omega * stage_elements(is_inductor_parallel_local).* stage_elementsQ(is_inductor_parallel_local) ));   
            if ~isempty(E_cp)
                YPm_Tx = YPm_Tx + diag(E_cp.*Cap_loss./(E_cp+Cap_loss));
                YPm_loss_Tx = YPm_loss_Tx + diag(E_cp);
            end
            if ~isempty(E_lp)
                YPm_Tx = YPm_Tx + diag(E_lp.*Ind_loss./(E_lp+Ind_loss));
                YPm_loss_Tx = YPm_loss_Tx + diag(E_lp);
            end
            if ~isempty(E_rp)
                YPm_Tx = YPm_Tx + E_rp;
            end
        end
        if sum(is_capacitor_series_local)> 0 || sum(is_inductor_series_local)>0 || sum(is_resistor_series_local)>0
            E_cs = (1./(1i*omega*stage_elements(is_capacitor_series_local))); 
            E_ls = (1i*omega*stage_elements(is_inductor_series_local));    
            E_rs = diag(stage_elements(is_resistor_series_local));    
            
            Cap_loss  = (1./(omega * stage_elements(is_capacitor_series_local).* stage_elementsQ(is_capacitor_series_local) ));   
            Ind_loss  = (omega * stage_elements(is_inductor_series_local) .* stage_elementsQ(is_inductor_series_local) ); 
            ZPm_Tx = inv(YPm_Tx);
            ZPm_loss_Tx = inv(YPm_loss_Tx);
            if ~isempty(E_cs)
                ZPm_Tx = ZPm_Tx + diag(E_cs.*Cap_loss./(E_cs+Cap_loss));
                ZPm_loss_Tx = ZPm_loss_Tx + diag(E_cs);
            end
            if ~isempty(E_ls)
                ZPm_Tx = ZPm_Tx + diag(E_ls.*Ind_loss./(E_ls+Ind_loss));
                ZPm_loss_Tx = ZPm_loss_Tx + diag(E_ls);
            end
            if ~isempty(E_rs)
                ZPm_Tx = ZPm_Tx + E_rs;
            end
            YPm_Tx = inv(ZPm_Tx);
            YPm_loss_Tx = inv(ZPm_loss_Tx);
        end
    end
    ZPm = inv(YPm_Tx);
    SPm = np_z2s(inv(YPm_Tx),z0);

    M_cal.rx = (M_decouple_Rx)*M_match_Rx;
    M_cal.tx = M_detune_Tx*calibration_matching(SPs,SPm,z0);
    
    ZPm_loss = inv(YPm_loss_Tx);
    phi_lumped_elements.Tx = diag(1/2*diag(abs(real(ZPm_loss-ZPm))).*(1./abs(diag(ZPm))).^2);
    phi_lumped_elements.Rx = 1/2*diag(abs(real(ZP_check_loss-ZP_check))).*(1./abs(diag(ZP_check))).^2;

end

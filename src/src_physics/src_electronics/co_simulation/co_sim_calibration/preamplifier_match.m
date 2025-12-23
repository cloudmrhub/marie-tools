function [M_match_Rx,...
          M_decouple_Rx,...
          SP_check,...
          ZP_check,...
          ZP_check_loss] = preamplifier_match(YPm,...
                                              YPm_loss,...
                                              preamp_res,...
                                              z0,...
                                              length_matching_network,...
                                              port_order,...
                                              is_capacitor_parallel,...
                                              is_inductor_parallel,...
                                              is_resistor_parallel,...
                                              is_capacitor_series,...
                                              is_inductor_series,...
                                              is_resistor_series,...
                                              matching_elements,...
                                              Q_val,...
                                              omega)

    P             = size(YPm,1);
    SP_check      = zeros(P,1);
    ZP_check      = zeros(P,1);
    ZP_check_loss = zeros(P,1);
    M_match_Rx    = zeros(P,P);
    M_decouple_Rx = zeros(P,P);
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
        M_decouple_Rx(:, k) = decouple_col;
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
            stage_elementsQ(good) = Q_val(ii(good));
            is_capacitor_parallel_local = is_capacitor_parallel(port_order==i);
            is_inductor_parallel_local  = is_inductor_parallel(port_order==i);
            is_resistor_parallel_local  = is_resistor_parallel(port_order==i);
            is_capacitor_series_local   = is_capacitor_series(port_order==i);
            is_inductor_series_local    = is_inductor_series(port_order==i);
            is_resistor_series_local    = is_resistor_series(port_order==i);
            if sum(is_capacitor_parallel_local)> 0 || sum(is_inductor_parallel_local)>0 || sum(is_resistor_parallel_local)>0
                E_cp = diag(1./(1i* omega * stage_elements(is_capacitor_parallel_local))); 
                E_lp = diag(1i*omega*stage_elements(is_inductor_parallel_local));    
                E_rp = diag(stage_elements(is_resistor_parallel_local));   
            
                Cap_loss  = diag( 1./( omega * stage_elements(is_capacitor_parallel_local) .* stage_elementsQ(is_capacitor_parallel_local) ));
                Ind_loss  = diag( omega * stage_elements(is_inductor_parallel_local)./ stage_elementsQ(is_inductor_parallel_local) );    
                if ~isempty(E_cp)
                    YPm_temp = YPm_temp + 1./(E_cp(k,k)+Cap_loss(k,k));
                    YPm_tempQ = YPm_tempQ + 1./E_cp(k,k);
                end
                if ~isempty(E_lp)
                    YPm_temp = YPm_temp + 1./(E_lp(k,k)+Ind_loss(k,k));
                    YPm_tempQ = YPm_tempQ + 1./E_lp(k,k);
                end
                if ~isempty(E_rp)
                    YPm_temp = YPm_temp + 1./E_rp(k,k);
                end
            end
            if sum(is_capacitor_series_local)> 0 || sum(is_inductor_series_local)>0 || sum(is_resistor_series_local)>0
                E_cs = diag(1./(1i*omega*stage_elements(is_capacitor_series_local))); 
                E_ls = diag(1i*omega*stage_elements(is_inductor_series_local));    
                E_rs = diag(stage_elements(is_resistor_series_local));    
            
                Cap_loss  = diag( 1./(omega * stage_elements(is_capacitor_series_local).* stage_elementsQ(is_capacitor_series_local) ));   
                Ind_loss  = diag( omega * stage_elements(is_inductor_series_local) ./ stage_elementsQ(is_inductor_series_local) ); 
                ZPm_temp = inv(YPm_temp);
                ZPm_tempQ = inv(YPm_tempQ);
                if ~isempty(E_cs)
                    ZPm_temp = ZPm_temp + E_cs(k,k)+Cap_loss(k,k);
                    ZPm_tempQ = ZPm_tempQ + E_cs(k,k);
                end
                if ~isempty(E_ls)
                    ZPm_temp = ZPm_temp + E_ls(k,k)+Ind_loss(k,k);
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

end
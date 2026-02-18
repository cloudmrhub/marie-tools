function [M_match_Rx,...
          M_decouple_Rx,...
          SP_check,...
          ZP_check,...
          ZP_check_loss] = preamplifier_match(YPm,...
                                              YPm_loss,...
                                              preamp_res,...
                                              z0,...
                                              ~,...
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
    
    % Get Matching Slots
    patched = zeros(size(port_order));   % preallocate
    patched(port_order ~= 0) = matching_elements;
    matching_elements = patched;
    patched = zeros(size(port_order));   % preallocate
    patched(port_order ~= 0) = Q_val;
    Q_val = patched;

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

        % Matching Elements
        idx = find(port_order == 1);
        start_idx = idx(k);
        if length(idx) >= k+1
            end_idx = idx(k+1) - 1;
        else
            end_idx = numel(port_order);
        end
        segment                     = nonzeros(port_order(start_idx:end_idx));
        local_matching_elements     = nonzeros(matching_elements(start_idx:end_idx));
        stage_elementsQ             = nonzeros(Q_val(start_idx:end_idx));
        is_capacitor_parallel_local = is_capacitor_parallel(start_idx:end_idx);
        is_inductor_parallel_local  = is_inductor_parallel(start_idx:end_idx);
        is_resistor_parallel_local  = is_resistor_parallel(start_idx:end_idx);
        is_capacitor_series_local   = is_capacitor_series(start_idx:end_idx);
        is_inductor_series_local    = is_inductor_series(start_idx:end_idx);
        is_resistor_series_local    = is_resistor_series(start_idx:end_idx);

        for i = 1:length(segment)
            icpl = is_capacitor_parallel_local(i);
            iipl = is_inductor_parallel_local(i); 
            irpl = is_resistor_parallel_local(i); 
            icsl = is_capacitor_series_local(i); 
            iisl = is_inductor_series_local(i); 
            irsl = is_resistor_series_local(i);   
            if icpl 
                E_cp = 1/(1i* omega * local_matching_elements(i)); 
                Cap_loss  = 1/( omega * local_matching_elements(i) * stage_elementsQ(i));   
                YPm_temp = YPm_temp + 1/(E_cp+Cap_loss);
                YPm_tempQ = YPm_tempQ + 1/E_cp;
            elseif iipl
                E_lp = 1i*omega*local_matching_elements(i);    
                Ind_loss  = omega * local_matching_elements(i)./ stage_elementsQ(i);    
                YPm_temp = YPm_temp + 1/(E_lp+Ind_loss);
                YPm_tempQ = YPm_tempQ + 1/E_lp;
            elseif irpl 
                E_rp = local_matching_elements(i);   
                YPm_temp = YPm_temp + 1/E_rp;
            elseif icsl
                ZPm_temp = inv(YPm_temp); 
                ZPm_tempQ = inv(YPm_tempQ);
                E_cs = 1/(1i*omega*local_matching_elements(i)); 
                Cap_loss  = 1/(omega * local_matching_elements(i) * stage_elementsQ(i));   
                ZPm_temp = ZPm_temp + E_cs+Cap_loss;
                ZPm_tempQ = ZPm_tempQ + E_cs;
                YPm_temp = inv(ZPm_temp);
                YPm_tempQ = inv(ZPm_tempQ);
            elseif iisl
                ZPm_temp = inv(YPm_temp); 
                ZPm_tempQ = inv(YPm_tempQ);
                E_ls = 1i*omega*local_matching_elements(i);    
                Ind_loss  = omega * local_matching_elements(i) / stage_elementsQ(i); 
                ZPm_temp = ZPm_temp + E_ls+Ind_loss;
                ZPm_tempQ = ZPm_tempQ + E_ls;
                YPm_temp = inv(ZPm_temp);
                YPm_tempQ = inv(ZPm_tempQ);
            elseif irsl
                ZPm_temp = inv(YPm_temp); 
                ZPm_tempQ = inv(YPm_tempQ);
                E_rs = local_matching_elements(i);    
                ZPm_temp = ZPm_temp + E_rs;
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
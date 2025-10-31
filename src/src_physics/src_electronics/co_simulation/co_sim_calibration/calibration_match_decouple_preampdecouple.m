function [M_cal,YPm,ZPm,SPm,SP_check] = calibration_match_decouple_preampdecouple(matching_elements,YPm,matching_load_string,port_order,omega,z0,preamp_res)
    
    is_capacitor_parallel = strcmpi(matching_load_string, 'capacitorParallel');
    is_inductor_parallel  = strcmpi(matching_load_string, 'inductorParallel');
    is_resistor_parallel  = strcmpi(matching_load_string, 'resistorParallel');
    is_capacitor_series   = strcmpi(matching_load_string, 'capacitorSeries');
    is_inductor_series    = strcmpi(matching_load_string, 'inductorSeries');
    is_resistor_series    = strcmpi(matching_load_string, 'resistorSeries');

    length_matching_network = unique(max(port_order));

    % Model preamplifier effect
    P = size(YPm,1);
    SP_check = zeros(P,1);
    M_match_Rx = zeros(P,P);
    M_decouple_Rx = zeros(P,P);
    for k = 1:P
        rps = setdiff(1:P, k);
        YPn_temp = YPm;
         % Decouple
        YPn_temp(rps, rps) = YPn_temp(rps, rps) + 1/preamp_res;
        ZPn_temp = inv(YPn_temp);
        decouple_col = zeros(P,1);
        decouple_col(k,1) = ZPn_temp(k,k);
        decouple_col(rps,1) = ZPn_temp(rps,k);
        decouple_col = decouple_col / ZPn_temp(k,k);
        M_decouple_Rx(:, k) = decouple_col;
        % Matching Elements
        YPm_temp = YPn_temp(k,k)-YPn_temp(k,rps)/YPn_temp(rps,rps)*YPn_temp(rps,k);
        SPs = np_z2s(inv(YPm_temp),z0);
        for i = 1:length_matching_network
            stage_mask           = (port_order == i);
            existsXx             = is_capacitor_parallel | is_inductor_parallel | is_resistor_parallel | is_capacitor_series | is_inductor_series | is_resistor_series;
            slot2idx             = zeros(numel(existsXx),1);
            slot2idx(existsXx)   = 1:nnz(existsXx);          
            ii                   = slot2idx(stage_mask); 
            stage_elements       = zeros(nnz(stage_mask),1);
            good                 = ii > 0;
            stage_elements(good) = matching_elements(ii(good));
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
                if ~isempty(E_cp)
                    YPm_temp = YPm_temp + E_cp(k,k);
                end
                if ~isempty(E_lp)
                    YPm_temp = YPm_temp + E_lp(k,k);
                end
                if ~isempty(E_rp)
                    YPm_temp = YPm_temp + E_rp(k,k);
                end
            end
            if sum(is_capacitor_series_local)> 0 || sum(is_inductor_series_local)>0 || sum(is_resistor_series_local)>0
                E_cs = diag(1./(1i*omega*stage_elements(is_capacitor_series_local))); 
                E_ls = diag(1i*omega*stage_elements(is_inductor_series_local));    
                E_rs = diag(stage_elements(is_resistor_series_local));    
                ZPm_temp = inv(YPm_temp);
                if ~isempty(E_cs)
                    ZPm_temp = ZPm_temp + E_cs(k,k);
                end
                if ~isempty(E_ls)
                    ZPm_temp = ZPm_temp + E_ls(k,k);
                end
                if ~isempty(E_rs)
                    ZPm_temp = ZPm_temp + E_rs(k,k);
                end
                YPm_temp = inv(ZPm_temp);
            end
        end
        SPm = np_z2s(inv(YPm_temp),z0);
        M_match_Rx(k,k) = calibration_matching(SPs,SPm,z0);
        SP_check(k,1) = SPm;
    end

    SPs = np_z2s(inv(YPm),z0);
    % Model matching elements
    for i = 1:length_matching_network
        stage_mask           = (port_order == i);
        existsXx             = is_capacitor_parallel | is_inductor_parallel | is_resistor_parallel | is_capacitor_series | is_inductor_series | is_resistor_series;
        slot2idx             = zeros(numel(existsXx),1);
        slot2idx(existsXx)   = 1:nnz(existsXx);          
        ii                   = slot2idx(stage_mask); 
        stage_elements       = zeros(nnz(stage_mask),1);
        good                 = ii > 0;
        stage_elements(good) = matching_elements(ii(good));
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
            if ~isempty(E_cp)
                YPm = YPm + E_cp;
            end
            if ~isempty(E_lp)
                YPm = YPm + E_lp;
            end
            if ~isempty(E_rp)
                YPm = YPm + E_rp;
            end
        end
        if sum(is_capacitor_series_local)> 0 || sum(is_inductor_series_local)>0 || sum(is_resistor_series_local)>0
            E_cs = diag(1./(1i*omega*stage_elements(is_capacitor_series_local))); 
            E_ls = diag(1i*omega*stage_elements(is_inductor_series_local));    
            E_rs = diag(stage_elements(is_resistor_series_local));    
            ZPm = inv(YPm);
            if ~isempty(E_cs)
                ZPm = ZPm + E_cs;
            end
            if ~isempty(E_ls)
                ZPm = ZPm + E_ls;
            end
            if ~isempty(E_rs)
                ZPm = ZPm + E_rs;
            end
            YPm = inv(ZPm);
        end
    end
    ZPm = inv(YPm);
    SPm = np_z2s(inv(YPm),z0);

    M_cal.rx = M_decouple_Rx*M_match_Rx;
    M_cal.tx = calibration_matching(SPs,SPm,z0);

end

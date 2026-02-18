function [ZP_check] = preamplifier_match_optim(YPm,...
                                               preamp_res,...
                                               ~,...
                                               port_order,...
                                               is_capacitor_parallel,...
                                               is_inductor_parallel,...
                                               is_resistor_parallel,...
                                               is_capacitor_series,...
                                               is_inductor_series,...
                                               is_resistor_series,...
                                               matching_elements,...
                                               omega)

    % Model preamplifier effect
    P = size(YPm,1);
    ZP_check = zeros(P,1);
    
    % Get Matching Slots
    patched = zeros(size(port_order));   % preallocate
    patched(port_order ~= 0) = matching_elements;
    matching_elements = patched;

    for k = 1:P
        rps = setdiff(1:P,k);
        YPm_temp = YPm(k,k) - YPm(k,rps)/(YPm(rps,rps) + 1/preamp_res)*YPm(rps,k);

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
                E_cp = 1i* omega * local_matching_elements(i); 
                YPm_temp = YPm_temp + E_cp; 
            elseif iipl
                E_lp = 1./(1i*omega*local_matching_elements(i)); 
                YPm_temp = YPm_temp + E_lp;
            elseif irpl 
                E_rp = 1./(local_matching_elements(i)); 
                YPm_temp = YPm_temp + E_rp; 
            elseif icsl
                ZPm_temp = inv(YPm_temp); 
                E_cs = 1./(1i*omega*local_matching_elements(i)); 
                ZPm_temp = ZPm_temp + E_cs; 
                YPm_temp = inv(ZPm_temp);
            elseif iisl
                ZPm_temp = inv(YPm_temp); 
                E_ls = 1i*omega*local_matching_elements(i);
                ZPm_temp = ZPm_temp + E_ls; 
                YPm_temp = inv(ZPm_temp);
            elseif irsl
                ZPm_temp = inv(YPm_temp); 
                E_rs = local_matching_elements(i); 
                ZPm_temp = ZPm_temp + E_rs;
                YPm_temp = inv(ZPm_temp);
            end 
        end
        ZP_check(k,1) = inv(YPm_temp);
    end

end
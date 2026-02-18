function [ZPm] = match_optim(YPm,...
                             length_matching_network,...
                             port_order,...
                             is_capacitor_parallel,...
                             is_inductor_parallel,...
                             is_resistor_parallel,...
                             is_capacitor_series,...
                             is_inductor_series,...
                             is_resistor_series,...
                             matching_elements,...
                             omega)    


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
                YPm(logical(is_capacitor_parallel_local),logical(is_capacitor_parallel_local)) = YPm(logical(is_capacitor_parallel_local),logical(is_capacitor_parallel_local)) + E_cp;
            end
            if ~isempty(E_lp)
                YPm(logical(is_inductor_parallel_local),logical(is_inductor_parallel_local)) = YPm(logical(is_inductor_parallel_local),logical(is_inductor_parallel_local)) + E_lp;
            end
            if ~isempty(E_rp)
                YPm(logical(is_resistor_parallel_local),logical(is_resistor_parallel_local)) = YPm(logical(is_resistor_parallel_local),logical(is_resistor_parallel_local)) + E_rp;
            end
        end

        if sum(is_capacitor_series_local)> 0 || sum(is_inductor_series_local)>0 || sum(is_resistor_series_local)>0
            E_cs = diag(1./(1i*omega*stage_elements(is_capacitor_series_local))); 
            E_ls = diag(1i*omega*stage_elements(is_inductor_series_local));    
            E_rs = diag(stage_elements(is_resistor_series_local));    
            ZPm = inv(YPm);
            if ~isempty(E_cs)
                ZPm(logical(is_capacitor_series_local),logical(is_capacitor_series_local)) = ZPm(logical(is_capacitor_series_local),logical(is_capacitor_series_local)) + E_cs;
            end
            if ~isempty(E_ls)
                ZPm(logical(is_inductor_series_local),logical(is_inductor_series_local)) = ZPm(logical(is_inductor_series_local),logical(is_inductor_series_local)) + E_ls;
            end
            if ~isempty(E_rs)
                ZPm(logical(is_resistor_series_local),logical(is_resistor_series_local)) = ZPm(logical(is_resistor_series_local),logical(is_resistor_series_local)) + E_rs;
            end
            YPm = inv(ZPm);
        end

    end
    ZPm = inv(YPm);

end
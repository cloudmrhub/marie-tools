function [SPm,ZPm,ZPm_loss] = match(YPm,...
                                    YPm_loss,...
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
            E_cp = (1./(1i* omega * stage_elements(is_capacitor_parallel_local))); 
            E_lp = (1i*omega*stage_elements(is_inductor_parallel_local));    
            E_rp = stage_elements(is_resistor_parallel_local);   
            
            Cap_loss  = (1./(omega * stage_elements(is_capacitor_parallel_local) .* stage_elementsQ(is_capacitor_parallel_local) ));
            Ind_loss  = ( omega * stage_elements(is_inductor_parallel_local)./ stage_elementsQ(is_inductor_parallel_local) );    
            if ~isempty(E_cp)
                YPm = YPm + diag(1./(E_cp+Cap_loss));
                YPm_loss = YPm_loss + diag(1./E_cp);
            end
            if ~isempty(E_lp)
                YPm = YPm + diag(1./(E_lp+Ind_loss));
                YPm_loss = YPm_loss + diag(1./E_lp);
            end
            if ~isempty(E_rp)
                YPm = YPm + diag(1./E_rp);
            end
        end
        if sum(is_capacitor_series_local)> 0 || sum(is_inductor_series_local)>0 || sum(is_resistor_series_local)>0
            E_cs = (1./(1i*omega*stage_elements(is_capacitor_series_local))); 
            E_ls = (1i*omega*stage_elements(is_inductor_series_local));    
            E_rs = stage_elements(is_resistor_series_local);   
            
            Cap_loss  = ( 1./(omega * stage_elements(is_capacitor_series_local).* stage_elementsQ(is_capacitor_series_local) ));   
            Ind_loss  = ( omega * stage_elements(is_inductor_series_local) ./ stage_elementsQ(is_inductor_series_local) );  
            ZPm = inv(YPm);
            ZPm_loss = inv(YPm_loss);
            if ~isempty(E_cs)
                ZPm = ZPm + diag(E_cs + Cap_loss);
                ZPm_loss = ZPm_loss + diag(E_cs);
            end
            if ~isempty(E_ls)
                ZPm = ZPm + diag(E_ls+Ind_loss);
                ZPm_loss = ZPm_loss + diag(E_ls);
            end
            if ~isempty(E_rs)
                ZPm = ZPm + diag(E_rs);
            end
            YPm = inv(ZPm);
            YPm_loss = inv(ZPm_loss);
        end
    end
    ZPm = inv(YPm);
    SPm = np_z2s(inv(YPm),z0);
    ZPm_loss = inv(YPm_loss);

end
function [M_cal,...
          phi_lumped_elements,...
          YPm,...
          ZPm,...
          SPm] = calibration_match_decouple(matching_elements,...
                                            YPm,...
                                            matching_load_string,...
                                            port_order,...
                                            omega,...
                                            z0,...
                                            Q_val)
    
    is_capacitor_parallel = strcmpi(matching_load_string, 'capacitorParallel');
    is_inductor_parallel  = strcmpi(matching_load_string, 'inductorParallel');
    is_resistor_parallel  = strcmpi(matching_load_string, 'resistorParallel');
    is_capacitor_series   = strcmpi(matching_load_string, 'capacitorSeries');
    is_inductor_series    = strcmpi(matching_load_string, 'inductorSeries');
    is_resistor_series    = strcmpi(matching_load_string, 'resistorSeries');

    length_matching_network = unique(max(port_order));
    SPs = np_z2s(inv(YPm),z0);

    YPm_loss = YPm;
    % Model matching elements
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
            E_cp = (1i* omega * stage_elements(is_capacitor_parallel_local)); 
            E_lp = (1./(1i*omega*stage_elements(is_inductor_parallel_local)));    
            E_rp = diag(1./(stage_elements(is_resistor_parallel_local)));   
            
            Cap_loss  = (omega * stage_elements(is_capacitor_parallel_local) .* stage_elementsQ(is_capacitor_parallel_local) );
            Ind_loss  = (1./(omega * stage_elements(is_inductor_parallel_local).* stage_elementsQ(is_inductor_parallel_local) ));    
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
            
            Cap_loss  = (1./(omega * stage_elements(is_capacitor_series_local).* stage_elementsQ(is_capacitor_series_local) ));   
            Ind_loss  = (omega * stage_elements(is_inductor_series_local) .* stage_elementsQ(is_inductor_series_local) ); 
            ZPm = inv(YPm);
            ZPm_loss = inv(YPm_loss);
            if ~isempty(E_cs)
                ZPm = ZPm + diag(E_cs.*Cap_loss./(E_cs+Cap_loss));
                ZPm_loss = ZPm_loss + E_cs;
            end
            if ~isempty(E_ls)
                ZPm = ZPm + diag(E_ls.*Ind_loss./(E_ls+Ind_loss));
                ZPm_loss = ZPm_loss + E_ls;
            end
            if ~isempty(E_rs)
                ZPm = ZPm + E_rs;
            end
            YPm = inv(ZPm);
            YPm_loss = inv(ZPm_loss);
        end
    end
    ZPm = inv(YPm);
    SPm = np_z2s(inv(YPm),z0);
    M_cal = calibration_matching(SPs,SPm,z0);

    ZPm_loss = inv(YPm_loss);
    phi_lumped_elements.Tx = 1/2*diag(abs(real(ZPm_loss-ZPm))).*(1./abs(diag(ZPm))).^2;
    phi_lumped_elements.Rx = 1/2*diag(abs(real(ZPm_loss-ZPm))).*(1./abs(diag(ZPm))).^2;

end

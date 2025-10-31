function[figure_freq] = get_ZPm_TxRx_freq(X,...
                                          YPn,...
                                          symmetries,...
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
                                          z0, ...
                                          preamp_res)

    freq                 = omega/(2*pi);
    freq_range           = linspace(freq * 0.8, freq * 1.2, 5000); 
    [~, freq_idx]        = min(abs(freq_range - freq)); 
    freq_range(freq_idx) = freq; 
    omega_range          = freq_range*2*pi;
      
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
    length_matching_network            = unique(max(port_order));
    M_detailed                         = nonzeros(matching_mask_detailed);
    F_detailed                         = nonzeros(fixed_matching_mask_detailed);
    N_detailed                         = nonzeros(tuning_mask_detailed);
    [M_detailed,F_detailed,N_detailed] = map_to_local(M_detailed,F_detailed,N_detailed);
    MN_detailed                        = sort([M_detailed(:); N_detailed(:)]);
    MF_detailed                        = sort([M_detailed(:); F_detailed(:)]);

    lumped_elements                  = X(symmetries); 
    all_elements                     = zeros(numel(M_detailed)+numel(F_detailed)+numel(N_detailed),1);
    all_elements(MN_detailed)        = lumped_elements;
    all_elements(F_detailed)         = fixed_matching_elements(F_detailed);
    matching_elements                = all_elements(MF_detailed);
    true_tuning_elements             = all_elements(N_detailed);
    is_capacitor_parallel            = strcmpi(matching_load_string, 'capacitorParallel');
    is_inductor_parallel             = strcmpi(matching_load_string, 'inductorParallel');
    is_resistor_parallel             = strcmpi(matching_load_string, 'resistorParallel');
    is_capacitor_series              = strcmpi(matching_load_string, 'capacitorSeries');
    is_inductor_series               = strcmpi(matching_load_string, 'inductorSeries');
    is_resistor_series               = strcmpi(matching_load_string, 'resistorSeries');
    is_capacitor                     = strcmpi(tuning_load_string, 'capacitor');
    is_inductor                      = strcmpi(tuning_load_string, 'inductor');
    is_resistor                      = strcmpi(tuning_load_string, 'resistor');
    is_mutual_inductor               = strcmpi(tuning_load_string, 'mutual_inductor');
    is_mutual_inductance_coefficient = strcmpi(tuning_load_string, 'mutual_inductance_coefficient');
    mutual_inductance_coefficients   = lumped_elements(end-sum(is_mutual_inductance_coefficient)+1:end);
    mutual_inductances               = true_tuning_elements(is_mutual_inductor);

    ZPm_plot_Rx = zeros(size(M,1),size(M,1),length(freq_range));
    SPm_plot_Rx = zeros(size(M,1),size(M,1),length(freq_range));
    ZPm_plot_Tx = zeros(size(M,1),size(M,1),length(freq_range));
    SPm_plot_Tx = zeros(size(M,1),size(M,1),length(freq_range));

    for freq_id = 1:length(freq_range)

        YPn_temp = YPn;

        % Tuning Elements
        Y_vals               = zeros(size(true_tuning_elements));
        Y_vals(is_capacitor) = 1i * omega_range(freq_id) .* true_tuning_elements(is_capacitor);
        Y_vals(is_inductor)  = 1 ./ (1i * omega_range(freq_id) .* true_tuning_elements(is_inductor));
        Y_vals(is_resistor)  = 1 ./ true_tuning_elements(is_resistor);
        Y_vals(is_mutual_inductor) = 1 ./ (1i*omega_range(freq_id) .* true_tuning_elements(is_mutual_inductor));
        
        % Construct diagonal admittance matrix and add values
        E_tu          = diag(Y_vals(:));
        YPn_temp(N,N) = YPn_temp(N,N) + E_tu;
        idx           = sub2ind(size(YPn_temp), D1, D2);
        YPn_temp(idx) = YPn_temp(idx) + 1 ./ (1i*omega_range(freq_id) .* mutual_inductance_coefficients.*sqrt(mutual_inductances(mapped_D1).*mutual_inductances(mapped_D2)));
        idx           = sub2ind(size(YPn_temp), D2, D1);
        YPn_temp(idx) = YPn_temp(idx) + 1 ./ (1i*omega_range(freq_id) .* mutual_inductance_coefficients.*sqrt(mutual_inductances(mapped_D2).*mutual_inductances(mapped_D1)));
        
        % Add tuning elements
        YPm = YPn_temp(M,M)-YPn_temp(M,N)/YPn_temp(N,N)*YPn_temp(N,M);
    
        % Model preamplifier effect
        P = size(YPm,1);
        ZP_check = zeros(P,1);
        SP_check = zeros(P,1);
        for k = 1:P
            rps = setdiff(1:P,k);
            YPm_temp = YPm;
            YPm_temp(rps,rps) = YPm_temp(rps,rps) + 1/preamp_res;
            % Matching Elements
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
                    E_cp = diag(1i* omega_range(freq_id) * stage_elements(is_capacitor_parallel_local)); 
                    E_lp = diag(1./(1i*omega_range(freq_id)*stage_elements(is_inductor_parallel_local)));    
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
                    E_cs = diag(1./(1i*omega_range(freq_id)*stage_elements(is_capacitor_series_local))); 
                    E_ls = diag(1i*omega_range(freq_id)*stage_elements(is_inductor_series_local));    
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
            SP_check(k,1) = np_z2s(ZP_check(k,1),z0);
        end
        ZPm_plot_Rx(:,:,freq_id) = diag(ZP_check);
        SPm_plot_Rx(:,:,freq_id) = diag(SP_check);

        % Matching Elements without Preamps
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
                E_cp = diag(1i* omega_range(freq_id) * stage_elements(is_capacitor_parallel_local)); 
                E_lp = diag(1./(1i*omega_range(freq_id)*stage_elements(is_inductor_parallel_local)));    
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
                E_cs = diag(1./(1i*omega_range(freq_id)*stage_elements(is_capacitor_series_local))); 
                E_ls = diag(1i*omega_range(freq_id)*stage_elements(is_inductor_series_local));    
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
        ZPm_plot_Tx(:,:,freq_id) = ZPm;
        SPm_plot_Tx(:,:,freq_id) = SPm;

    end

    figure_freq.ZPm_plot_Rx = ZPm_plot_Rx;
    figure_freq.SPm_plot_Rx = SPm_plot_Rx;
    figure_freq.ZPm_plot_Tx = ZPm_plot_Tx;
    figure_freq.SPm_plot_Tx = SPm_plot_Tx;
    figure_freq.freq_range  = freq_range;
    figure_freq.freq_idx    = freq_idx;

end
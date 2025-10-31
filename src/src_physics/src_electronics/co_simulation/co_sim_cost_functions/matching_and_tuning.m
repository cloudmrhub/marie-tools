function [cf] = matching_and_tuning(X,...
                                    YPn,...
                                    symmetries,...
                                    M_detailed,...
                                    F_detailed,...
                                    fixed_matching_elements,...
                                    N_detailed,...
                                    MN_detailed,...
                                    MF_detailed,...
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
                                    omega,...
                                    z0)

    lumped_elements                  = X(symmetries); 
    all_elements                     = zeros(numel(M_detailed)+numel(F_detailed)+numel(N_detailed),1);
    all_elements(MN_detailed)        = lumped_elements;
    all_elements(F_detailed)         = fixed_matching_elements;
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

    % Matching Elements
    for i = 1:length_matching_network

        stage_mask     = (port_order == i);
        stage_elements = matching_elements(stage_mask);
    
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

    % Final inversions
    ZPm = inv(YPm);
    cf  = norm(diag(real(ZPm(Final_M,Final_M))-z0),'fro') + norm(diag(imag(ZPm)),'fro');

    % if size(ZPm,1) > 1
    %     cf = cf + norm(tril(ZPm,-1));
    % end

end
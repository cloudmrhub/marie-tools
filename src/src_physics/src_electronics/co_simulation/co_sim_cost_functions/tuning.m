function [cf,ZPm] = tuning(X,...
                           YPn,...
                           symmetries,...
                           M,...
                           N,...
                           idx12,...
                           idx21,...
                           mapped_D1,...
                           mapped_D2,...
                           tuning_load_string,...
                           omega)

    tuning_elements                  = X(symmetries); 
    is_capacitor                     = strcmpi(tuning_load_string,'capacitor');
    is_inductor                      = strcmpi(tuning_load_string,'inductor');
    is_resistor                      = strcmpi(tuning_load_string,'resistor');
    is_mutual_inductor               = strcmpi(tuning_load_string,'mutual_inductor');
    is_mutual_inductance_coefficient = strcmpi(tuning_load_string,'mutual_inductance_coefficient');
    true_tuning_elements             = tuning_elements(1:end-sum(is_mutual_inductance_coefficient));
    mutual_inductance_coefficients   = tuning_elements(end-sum(is_mutual_inductance_coefficient)+1:end);
    mutual_inductances               = tuning_elements(is_mutual_inductor);

    % Preallocate admittance values
    Y_vals                     = zeros(size(true_tuning_elements));
    Y_vals(is_capacitor)       = 1i * omega .* true_tuning_elements(is_capacitor);
    Y_vals(is_inductor)        = 1 ./ (1i*omega .* true_tuning_elements(is_inductor));
    Y_vals(is_resistor)        = 1 ./ true_tuning_elements(is_resistor);
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
    ZPm = inv(YPm);

    cf = norm(diag(imag(ZPm)),'fro');

    % if size(ZPm,1) > 1
    %     cf = cf + norm(tril(ZPm,-1));
    % end

end
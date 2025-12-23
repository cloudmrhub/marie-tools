function [YPm,YPm_loss,M_tune] = tune(true_tuning_elements,...
                                      omega,...
                                      Q_true_tuning_elements,...
                                      is_capacitor,...
                                      is_inductor,...
                                      is_mutual_inductor,...
                                      is_resistor,...
                                      YPn,...
                                      D1,...
                                      D2,...
                                      mutual_inductance_coefficients,...
                                      mutual_inductances,...
                                      mapped_D1,...
                                      mapped_D2,...
                                      N,...
                                      M)   

    Y_vals               = zeros(size(true_tuning_elements));
    Y_vals_loss          = zeros(size(true_tuning_elements));

    Cap_loss             = 1./ (omega .* true_tuning_elements(is_capacitor) .* Q_true_tuning_elements(is_capacitor));
    Ind_loss             = omega .* true_tuning_elements(is_inductor) ./ Q_true_tuning_elements(is_inductor);
    MutInd_loss          = omega .* true_tuning_elements(is_mutual_inductor) ./ Q_true_tuning_elements(is_mutual_inductor);
    Cap_Value            = 1 ./ (1i * omega .* true_tuning_elements(is_capacitor));
    Ind_Value            = 1i * omega .* true_tuning_elements(is_inductor);
    MutInd_Value         = 1i*omega .* true_tuning_elements(is_mutual_inductor);
    
    Y_vals(is_capacitor)            = 1./(Cap_Value + Cap_loss);
    Y_vals(is_inductor)             = 1./(Ind_Value + Ind_loss);
    Y_vals(is_resistor)             = 1./true_tuning_elements(is_resistor);
    Y_vals(is_mutual_inductor)      = 1./(MutInd_Value + MutInd_loss);
    Y_vals_loss(is_capacitor)       = 1./Cap_Value;
    Y_vals_loss(is_inductor)        = 1./Ind_Value;
    Y_vals_loss(is_mutual_inductor) = 1./MutInd_Value;

    YPn_loss = YPn;
    
    % Construct diagonal admittance matrix and add values
    E_tu     = diag(Y_vals(:));
    YPn(N,N) = YPn(N,N) + E_tu;
    idx      = sub2ind(size(YPn), D1, D2);
    YPn(idx) = YPn(idx) + 1 ./ (1i*omega .* mutual_inductance_coefficients.*sqrt(mutual_inductances(mapped_D1).*mutual_inductances(mapped_D2)));
    idx      = sub2ind(size(YPn), D2, D1);
    YPn(idx) = YPn(idx) + 1 ./ (1i*omega .* mutual_inductance_coefficients.*sqrt(mutual_inductances(mapped_D2).*mutual_inductances(mapped_D1)));
    ZPn      = inv(YPn);

    % Construct diagonal loss admittance matrix and add values
    E_tu_loss     = diag(Y_vals_loss(:));
    YPn_loss(N,N) = YPn_loss(N,N) + E_tu_loss;
    idx           = sub2ind(size(YPn_loss), D1, D2);
    YPn_loss(idx) = YPn_loss(idx) + 1 ./ (1i*omega .* mutual_inductance_coefficients.*sqrt(mutual_inductances(mapped_D1).*mutual_inductances(mapped_D2)));
    idx           = sub2ind(size(YPn_loss), D2, D1);
    YPn_loss(idx) = YPn_loss(idx) + 1 ./ (1i*omega .* mutual_inductance_coefficients.*sqrt(mutual_inductances(mapped_D2).*mutual_inductances(mapped_D1)));

    M_tune                = zeros(length(N)+length(M),length(M));
    M_tune(M,1:length(M)) = ZPn(M,M);
    M_tune(N,1:length(M)) = ZPn(N,M);
    M_tune                = M_tune/ZPn(M,M);

    % Add tuning elements
    YPm = YPn(M,M)-YPn(M,N)/YPn(N,N)*YPn(N,M);

    % Add tuning elements to loss
    YPm_loss = YPn_loss(M,M)-YPn_loss(M,N)/YPn_loss(N,N)*YPn_loss(N,M);

end
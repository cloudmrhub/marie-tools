function[MREDM] = co_simulation(MREDM)
    
    if isempty(MREDM.WIE.coil) && isempty(MREDM.SIE.coil) && isempty(MREDM.SIE.shield)
        return
    end

    % Get constants
    emc = MREDM.emc;

    if MREDM.inputs.tmd == 1
   
        % Get Coil Structs
        [RLC_org,RLC_tmd,RLC_org_h,RLC_tmd_h,RLC_org_w,RLC_tmd_w,RLC_org_s,RLC_tmd_s,separator_counter,shield_file,wire_file,coil_file] = co_simulation_coil_optimize_settings(MREDM);

        % Get Y parameters
        YPn = MREDM.fields.netp.YP_s;

        % PSO Options
        [options_PSO_1, options_PSO_2] = co_simulation_optimization_settings();

        coil_types                   = arrayfun(@(x) x.excitation.TxRx, RLC_tmd, 'UniformOutput', false);
        coil_types                   = string(coil_types);
        coil_type_Rx                 = all(coil_types == "Rx");
        coil_type_Tx                 = all(coil_types == "Tx");
        coil_type_Tx_and_Rx          = all(ismember(coil_types, ["Rx", "Tx"])) && all(ismember(["Rx", "Tx"], coil_types));
        coil_type_TxRx               = all(coil_types == "TxRx");
        coil_type_TxRx_and_Rx        = all(ismember(coil_types, ["Rx", "TxRx"])) && all(ismember(["Rx", "TxRx"], coil_types));
        coil_type_TxRx_and_Tx        = all(ismember(coil_types, ["TxRx", "Tx"])) && all(ismember(["TxRx", "Tx"], coil_types));
        coil_type_TxRx_and_Tx_and_Rx = all(ismember(coil_types, ["TxRx", "Rx", "Tx"])) && all(ismember(["TxRx", "Rx", "Tx"], coil_types));

        % Get Losses
        Q_val = co_simulation_get_Q_values(RLC_tmd);

        % Decision Tree
        if coil_type_Rx
            [return_elements,phi_lumped_elements,M_cal,YPm,ZPm,SPm,SP_det,figure_freq] = co_simulation_Rx_optim(YPn,emc,RLC_org,RLC_tmd,Q_val,options_PSO_1,options_PSO_2);
        elseif coil_type_Tx
            [return_elements,phi_lumped_elements,M_cal,YPm,ZPm,SPm,SP_det,figure_freq] = co_simulation_Tx_optim(YPn,emc,RLC_org,RLC_tmd,Q_val,options_PSO_1,options_PSO_2);
        elseif coil_type_Tx_and_Rx
            [return_elements,phi_lumped_elements,M_cal,YPm,ZPm,SPm,SP_det,figure_freq] = co_simulation_Rx_and_Tx_optim(YPn,emc,RLC_org,RLC_tmd,Q_val,options_PSO_1,options_PSO_2);
        elseif coil_type_TxRx
            [return_elements,phi_lumped_elements,M_cal,YPm,ZPm,SPm,SP_det,figure_freq] = co_simulation_TxRx_optim(YPn,emc,RLC_org,RLC_tmd,Q_val,options_PSO_1,options_PSO_2);
        elseif coil_type_TxRx_and_Rx
            [return_elements,phi_lumped_elements,M_cal,YPm,ZPm,SPm,SP_det,figure_freq] = co_simulation_TxRx_and_Rx_optim(YPn,emc,RLC_org,RLC_tmd,Q_val,options_PSO_1,options_PSO_2);
        elseif coil_type_TxRx_and_Tx
            [return_elements,phi_lumped_elements,M_cal,YPm,ZPm,SPm,SP_det,figure_freq] = co_simulation_TxRx_and_Tx_optim(YPn,emc,RLC_org,RLC_tmd,Q_val,options_PSO_1,options_PSO_2);
        elseif coil_type_TxRx_and_Tx_and_Rx
            [return_elements,phi_lumped_elements,M_cal,YPm,ZPm,SPm,SP_det,figure_freq] = co_simulation_TxRx_and_Rx_and_Tx_optim(YPn,emc,RLC_org,RLC_tmd,Q_val,options_PSO_1,options_PSO_2);
        end
    
        % Store Updated elements
        MREDM = co_simulation_store_elements(MREDM,separator_counter,shield_file,wire_file,coil_file,RLC_org_h,RLC_tmd_h,RLC_org_w,RLC_tmd_w,RLC_org_s,RLC_tmd_s,return_elements);
        MREDM.fields.netp.figure_freq = figure_freq;

    elseif MREDM.inputs.tmd == 0
        % Find Coil
        [RLC_org,~,~,~,~,~] = co_simulation_coil_settings(MREDM);

        coil_types                   = arrayfun(@(x) x.excitation.TxRx, RLC_org, 'UniformOutput', false);
        coil_types                   = string(coil_types);
        coil_type_Rx                 = all(coil_types == "Rx");
        coil_type_Tx                 = all(coil_types == "Tx");
        coil_type_Tx_and_Rx          = all(ismember(coil_types, ["Rx", "Tx"])) && all(ismember(["Rx", "Tx"], coil_types));
        coil_type_TxRx               = all(coil_types == "TxRx");
        coil_type_TxRx_and_Rx        = all(ismember(coil_types, ["Rx", "TxRx"])) && all(ismember(["Rx", "TxRx"], coil_types));
        coil_type_TxRx_and_Tx        = all(ismember(coil_types, ["TxRx", "Tx"])) && all(ismember(["TxRx", "Tx"], coil_types));
        coil_type_TxRx_and_Tx_and_Rx = all(ismember(coil_types, ["TxRx", "Rx", "Tx"])) && all(ismember(["TxRx", "Rx", "Tx"], coil_types));

        % Get Y parameters
        YPn = MREDM.fields.netp.YP_s;

        % Get Losses
        Q_val = co_simulation_get_Q_values(RLC_org);

        if coil_type_Rx
            [M_cal,phi_lumped_elements,YPm,ZPm,SPm,SP_det] = co_simulation_Rx_no_optim(YPn,emc,RLC_org,Q_val);
        elseif coil_type_Tx
            [M_cal,phi_lumped_elements,YPm,ZPm,SPm,SP_det] = co_simulation_Tx_no_optim(YPn,emc,RLC_org,Q_val);
        elseif coil_type_Tx_and_Rx
            [M_cal,phi_lumped_elements,YPm,ZPm,SPm,SP_det] = co_simulation_Rx_and_Tx_no_optim(YPn,emc,RLC_org,Q_val);
        elseif coil_type_TxRx
            [M_cal,phi_lumped_elements,YPm,ZPm,SPm,SP_det] = co_simulation_TxRx_no_optim(YPn,emc,RLC_org,Q_val);
        elseif coil_type_TxRx_and_Rx
            [M_cal,phi_lumped_elements,YPm,ZPm,SPm,SP_det] = co_simulation_Rx_and_TxRx_no_optim(YPn,emc,RLC_org,Q_val);
        elseif coil_type_TxRx_and_Tx
            [M_cal,phi_lumped_elements,YPm,ZPm,SPm,SP_det] = co_simulation_TxRx_and_Tx_no_optim(YPn,emc,RLC_org,Q_val);
        elseif coil_type_TxRx_and_Tx_and_Rx
            [M_cal,phi_lumped_elements,YPm,ZPm,SPm,SP_det] = co_simulation_TxRx_and_TxRx_no_optim(YPn,emc,RLC_org,Q_val);
        end
        
    end

    MREDM = co_simulation_fields_calibration(MREDM, ...
                                             M_cal, ...
                                             coil_type_Rx, ...
                                             coil_type_Tx, ...
                                             coil_type_Tx_and_Rx, ...
                                             coil_type_TxRx, ...
                                             coil_type_TxRx_and_Rx, ...
                                             coil_type_TxRx_and_Tx, ...
                                             coil_type_TxRx_and_Tx_and_Rx);
    
    MREDM.fields.netp.YP_s = YPm;
    MREDM.fields.netp.ZP_s = ZPm;
    MREDM.fields.netp.SP_s = SPm;
    MREDM.fields.netp.SP_p = SP_det;
    MREDM.fields.netp.phi_lumped_elements = phi_lumped_elements;

end
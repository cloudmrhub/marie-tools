function[MREDM] = co_simulation_fields_calibration(MREDM, ...
                                                   M_cal, ...
                                                   coil_type_Rx, ...
                                                   coil_type_Tx, ...
                                                   coil_type_Tx_and_Rx, ...
                                                   coil_type_TxRx, ...
                                                   coil_type_TxRx_and_Rx, ...
                                                   coil_type_TxRx_and_Tx, ...
                                                   coil_type_TxRx_and_Tx_and_Rx)

    % Store Updated fields
    if coil_type_Rx
        MREDM.fields.JcbRx = MREDM.fields.Jcb*M_cal;
        MREDM.fields.JcbTx = []; 
    elseif coil_type_Tx
        MREDM.fields.JcbRx = [];
        MREDM.fields.JcbTx = MREDM.fields.Jcb*M_cal; 
    elseif coil_type_Tx_and_Rx
        MREDM.fields.JcbRx = MREDM.fields.Jcb*M_cal.rx;
        MREDM.fields.JcbTx = MREDM.fields.Jcb*M_cal.tx; 
    elseif coil_type_TxRx
        MREDM.fields.JcbRx = MREDM.fields.Jcb*M_cal.rx;
        MREDM.fields.JcbTx = MREDM.fields.Jcb*M_cal.tx; 
    elseif coil_type_TxRx_and_Rx
        MREDM.fields.JcbRx = MREDM.fields.Jcb*M_cal.rx;
        MREDM.fields.JcbTx = MREDM.fields.Jcb*M_cal.tx; 
    elseif coil_type_TxRx_and_Tx
        MREDM.fields.JcbRx = MREDM.fields.Jcb*M_cal.rx;
        MREDM.fields.JcbTx = MREDM.fields.Jcb*M_cal.tx; 
    elseif coil_type_TxRx_and_Tx_and_Rx
        MREDM.fields.JcbRx = MREDM.fields.Jcb*M_cal.rx;
        MREDM.fields.JcbTx = MREDM.fields.Jcb*M_cal.tx; 
    end

    % Check for MRGF
    if isfield(MREDM.fields,'a')
        if coil_type_Rx
            MREDM.fields.aRx = MREDM.fields.a*M_cal;
            MREDM.fields.aTx = []; 
        elseif coil_type_Tx
            MREDM.fields.aRx = [];
            MREDM.fields.aTx = MREDM.fields.a*M_cal; 
        elseif coil_type_Tx_and_Rx
            MREDM.fields.aRx = MREDM.fields.a*M_cal.rx;
            MREDM.fields.aTx = MREDM.fields.a*M_cal.tx; 
        elseif coil_type_TxRx
            MREDM.fields.aRx = MREDM.fields.a*M_cal.rx;
            MREDM.fields.aTx = MREDM.fields.a*M_cal.tx; 
        elseif coil_type_TxRx_and_Rx
            MREDM.fields.aRx = MREDM.fields.a*M_cal.rx;
            MREDM.fields.aTx = MREDM.fields.a*M_cal.tx; 
        elseif coil_type_TxRx_and_Tx
            MREDM.fields.aRx = MREDM.fields.a*M_cal.rx;
            MREDM.fields.aTx = MREDM.fields.a*M_cal.tx; 
        elseif coil_type_TxRx_and_Tx_and_Rx
            MREDM.fields.aRx = MREDM.fields.a*M_cal.rx;
            MREDM.fields.aTx = MREDM.fields.a*M_cal.tx; 
        end
    end

end
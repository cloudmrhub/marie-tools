function[MREDM] = em_ehfield_wsvie(MREDM)

    try 
        gpu_flag = 0;
        MREDM.fields.e_field = MREDM.functions.e_field_wsvie(MREDM,gpu_flag);
        MREDM.fields.h_field = MREDM.functions.h_field_wsvie(MREDM,gpu_flag);
    catch
        gpu_flag = 5;
        warning('Out of GPU memory. Running in CPU.');
        MREDM.fields.e_field = MREDM.functions.e_field_wsvie(MREDM,gpu_flag);
        MREDM.fields.h_field = MREDM.functions.h_field_wsvie(MREDM,gpu_flag);
    end

    if ~(isempty(MREDM.WIE.coil) && isempty(MREDM.SIE.coil) && isempty(MREDM.SIE.shield))
        MREDM.fields.snr_field = em_SNR(MREDM);
        MREDM.fields.txe_field = em_TXE(MREDM);
    end
    
end
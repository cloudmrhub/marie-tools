function[MREDM] = em_ehfield_wsvie_mrgf(MREDM,Ue,Ub)

    MREDM.fields.e_field   = em_efield_svie_mrgf(MREDM,Ue);
    MREDM.fields.h_field   = em_hfield_svie_mrgf(MREDM,Ub);
    MREDM.fields.snr_field = em_SNR(MREDM);
    MREDM.fields.txe_field = em_TXE(MREDM);

end
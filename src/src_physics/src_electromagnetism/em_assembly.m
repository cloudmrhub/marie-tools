function[EP] = em_assembly(MREDM,geo_flag)

    idxS     = MREDM.dimensions.idxS;
    n1       = MREDM.dimensions.n1;
    n2       = MREDM.dimensions.n2;
    n3       = MREDM.dimensions.n3;
    Mc_inv   = zeros(n1,n2,n3);
    Mcr_inv  = zeros(n1,n2,n3);
    
    S  = load(MREDM.inputs.rhbm);
    er = S.RHBM.epsilon_r;
    se = S.RHBM.sigma_e;
    
    Mr            = er + se / (MREDM.emc.ce);
    Mc            = Mr - 1;
    Mcr           = Mc./ Mr;
    Mr_inv        = 1./Mr;
    Mc_inv(idxS)  = 1./Mc(idxS);
    Mcr_inv(idxS) = 1./Mcr(idxS);
    
    EP.er      = er;
    EP.se      = se;
    EP.Mr      = Mr;
    EP.Mc      = Mc;
    EP.Mcr     = Mcr;
    EP.Mr_inv  = Mr_inv;
    EP.Mc_inv  = Mc_inv;
    EP.Mcr_inv = Mcr_inv;

    if MREDM.inputs.pFFT_flag && geo_flag == 1 && (~isempty(MREDM.inputs.coil) || ~isempty(MREDM.inputs.wire))
        pfft_idxS  = MREDM.dimensions.pfft_idxS;
        pfft_n1    = MREDM.dimensions.pfft_n1;
        pfft_n2    = MREDM.dimensions.pfft_n2;
        pfft_n3    = MREDM.dimensions.pfft_n3;
        ext_Mc_inv = zeros(pfft_n1,pfft_n2,pfft_n3);
        
        ext_Mc_inv(pfft_idxS) = 1./(er(idxS) + se(idxS)/(MREDM.emc.ce) - 1);
        EP.ext_Mc_inv = ext_Mc_inv;
    end
       
end
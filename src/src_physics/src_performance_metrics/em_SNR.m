function [snr_field] = em_SNR(MREDM)
    dims            = MREDM.dimensions;
    n1              = dims.n1;
    n2              = dims.n2;
    n3              = dims.n3;
    idxS            = dims.idxS;
    N_scat          = dims.N_scat;
    res             = dims.res;
    omega           = MREDM.emc.omega;
    k_B             = MREDM.emc.k_B;
    T0              = MREDM.emc.T0;
    M_0             = MREDM.emc.M0;
    snr_num_scaling = res^3*omega*M_0/sqrt(4*k_B*T0);

    if isfield(MREDM.fields.h_field,'Rx')
        B1m_3D = reshape(MREDM.fields.h_field.Rx.b1m,n1*n2*n3,[]);
        phimat = MREDM.fields.e_field.Rx.phi;
        phimat_coil = MREDM.fields.e_field.Rx.phi_coil;
        phimat_lump = MREDM.fields.e_field.Rx.phi_lump;
        if isfield(MREDM.fields.netp,'phi_lumped_elements')
            phimat_lump = phimat_lump + MREDM.fields.netp.phi_lumped_elements.Rx;
        end
    else
        B1m_3D = reshape(MREDM.fields.h_field.Tx.b1m,n1*n2*n3,[]);
        phimat = MREDM.fields.e_field.Tx.phi;
        phimat_coil = MREDM.fields.e_field.Tx.phi_coil;
        phimat_lump = MREDM.fields.e_field.Tx.phi_lump;
        if isfield(MREDM.fields.netp,'phi_lumped_elements')
            phimat_lump = phimat_lump + MREDM.fields.netp.phi_lumped_elements.Tx;
        end
    end
    A_ports = size(B1m_3D,2);

    B1m = zeros(N_scat,A_ports);
    for port = 1:A_ports
        B1m_3D1 = B1m_3D(:,port);
        B1m(:,port) = B1m_3D1(idxS);
    end

    SNR = zeros(N_scat,1);
    invphimat = inv(phimat+phimat_coil+phimat_lump);
    for jj = 1:N_scat
        Svox         = transpose(B1m(jj,:));                
        SphiS        = Svox'*(invphimat*Svox); %#ok<MINV>
        invSphiS     = inv(SphiS);
        snr_map      = 1/sqrt(invSphiS);
        SNR(jj,1)    = snr_num_scaling * snr_map;
    end
    SNR_avg           = mean(nonzeros(SNR(:)));
    SNR_full          = zeros(n1,n2,n3);
    SNR_full(idxS)    = SNR;

    snr_field.SNR_avg = SNR_avg;
    snr_field.SNR     = SNR_full;

end
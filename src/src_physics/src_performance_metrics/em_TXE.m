function [txe_field] = em_TXE(MREDM)
    dims            = MREDM.dimensions;
    n1              = dims.n1;
    n2              = dims.n2;
    n3              = dims.n3;
    idxS            = dims.idxS;
    N_scat          = dims.N_scat;

    if isfield(MREDM.fields.h_field,'Tx')
        B1p_3D = reshape(MREDM.fields.h_field.Tx.b1p,n1*n2*n3,[]);
        phimat = MREDM.fields.e_field.Tx.phi;
        phimat_coil = MREDM.fields.e_field.Tx.phi_coil;
        phimat_lump = 0;%MREDM.fields.e_field.Tx.phi_lump;
    else
        B1p_3D = reshape(MREDM.fields.h_field.Rx.b1p,n1*n2*n3,[]);
        phimat = MREDM.fields.e_field.Rx.phi;
        phimat_coil = MREDM.fields.e_field.Rx.phi_coil;
        phimat_lump = 0;%MREDM.fields.e_field.Rx.phi_lump;
    end
    A_ports = size(B1p_3D,2);

    B1p = zeros(N_scat,A_ports);
    for port = 1:A_ports
        B1p_3D1 = B1p_3D(:,port);
        B1p(:,port) = B1p_3D1(idxS);
    end

    TXE = zeros(N_scat,1);
    invphimat = inv(phimat+phimat_coil+phimat_lump);
    for jj = 1:N_scat
        AA = conj(invphimat) * B1p(jj,:)'*B1p(jj,:);
        if ~isscalar(AA)
            TXE(jj,1) = eigs(AA,1);
        else
            TXE(jj,1) = AA;
        end
    end
    TXE_avg           = mean(nonzeros(TXE(:)));
    TXE_full          = zeros(n1,n2,n3);
    TXE_full(idxS)    = TXE;

    txe_field.TXE_avg = TXE_avg;
    txe_field.TXE     = TXE_full;

end
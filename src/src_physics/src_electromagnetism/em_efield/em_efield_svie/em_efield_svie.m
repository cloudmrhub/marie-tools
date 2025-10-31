function[e_field] = em_efield_svie(MREDM,gpu_flag)

    func         = MREDM.functions;
    dims         = MREDM.dimensions;
    ce           = MREDM.emc.ce;
    EP           = MREDM.EP;
    N_sie        = dims.N_sie;
    n1           = dims.n1;
    n2           = dims.n2;
    n3           = dims.n3;
    ql           = dims.ql;
    l            = dims.l;
    mask         = dims.mask;
    idx_solution = dims.idx_solution;
    res          = dims.res;
    idxS         = dims.idxS;
    max_rank_sie = dims.max_rank_N_sie;
    Mpow         = repmat(EP.Mc_inv(idxS),[ql 1])/ce; 
    fN           = to_GPU(MREDM.operators.N,gpu_flag);

    % Compute
    for i = 1:2
        if ~isempty(MREDM.fields.JcbTx) && i == 1
            Jb      = MREDM.fields.JcbTx(N_sie+1:end,:);
            Jc      = MREDM.fields.JcbTx(1:N_sie,:);
            A_ports = size(MREDM.fields.JcbTx,2);
        elseif ~isempty(MREDM.fields.JcbRx) && i == 2        
            Jb      = MREDM.fields.JcbRx(N_sie+1:end,:);
            Jc      = MREDM.fields.JcbRx(1:N_sie,:);
            A_ports = size(MREDM.fields.JcbRx,2);
        else
            continue
        end

        e_b          = zeros(n1,n2,n3,ql,A_ports);
        pabs         = zeros(n1,n2,n3,A_ports);
        sc_pabs      = zeros(A_ports,1);
        sc_pext      = zeros(A_ports,1);
        sc_psca      = zeros(A_ports,1);

        % Total Field and Power
        for port = 1:A_ports
            tested_einc           = func.mvp_Zbc(MREDM.operators.Zbc_N,Jc(:,port),n1,n2,n3,N_sie,max_rank_sie); 
            Jb_port               = zeros(n1,n2,n3,ql);
            Jb_port(idx_solution) = Jb(:,port);
            tested_esca           = gather(1/ce * (func.mvp_N(gpuArray(Jb_port),fN) - func.mvp_G(gpuArray(Jb_port),res)));
            eb_c                  = mask.*func.mvp_invG(reshape(squeeze(tested_einc),n1,n2,n3,ql),res);
            eb_b                  = mask.*func.mvp_invG(tested_esca,res);
            e_b(:,:,:,:,port)     = eb_c+eb_b;
            pabs(:,:,:,port)      = 1/2 * EP.se.* sqrt(abs(e_b(:,:,:,1,port)).^2+abs(e_b(:,:,:,1+l,port)).^2+abs(e_b(:,:,:,1+2*l,port)).^2).^2;
            G_J                   = func.mvp_G(gpuArray(Jb_port),res);
            Gx                    = gather(G_J(idx_solution));
            tested_einc           = ce*EP.Mc.*reshape(squeeze(tested_einc),n1,n2,n3,ql);
            tested_einc           = tested_einc(idx_solution);
            sc_pabs(port)         = 1/2 * real(Jb_port(idx_solution)'  * (Mpow .* Gx ));    
            sc_pext(port)         = -1/2 * real(Jb_port(idx_solution)' * (Mpow .* tested_einc(:)));
            sc_psca(port)         = -1/2*real((Jb_port(:))' * tested_esca(:));
        end
        
        % Body Noise Covariance Matrix
        Gx = zeros(size(Jb));
        for port = 1:A_ports 
            y               = zeros(n1,n2,n3,ql);
            y(idx_solution) = Jb(:,port);
            G_J             = func.mvp_G(gpuArray(y),res);
            Gx(:,port)      = gather(G_J(idx_solution));   
        end   
        phi = Jb' * (real(Mpow) .* Gx);
        phi = conj(phi);
    
        % Coil Noise Covariance Matrix
        phi_coil = zeros(A_ports,A_ports);
        phi_lump = zeros(A_ports,A_ports);
        for port = 1:A_ports 
            phi_coil(port,port) = 1/2 * Jc(:,port)' * (MREDM.SIE.Z_copper_loss * Jc(:,port));
            phi_lump(port,port) = 1/2 * Jc(:,port)' * (MREDM.SIE.Z_lumped_loss * Jc(:,port));
        end
        phi_coil = conj(phi_coil);
        phi_lump = conj(phi_lump);

        % Save to struct and return
        if ~isempty(MREDM.fields.JcbTx) && i == 1
            e_field.Tx.e_b      = e_b;
            e_field.Tx.pabs     = pabs;
            e_field.Tx.sc_pabs  = sc_pabs;
            e_field.Tx.sc_psca  = sc_psca;
            e_field.Tx.sc_pext  = sc_pext;
            e_field.Tx.phi      = phi;
            e_field.Tx.phi_coil = phi_coil;
            e_field.Tx.phi_lump = phi_lump;
        elseif ~isempty(MREDM.fields.JcbRx) && i == 2        
            e_field.Rx.e_b      = e_b;
            e_field.Rx.pabs     = pabs;
            e_field.Rx.sc_pabs  = sc_pabs;
            e_field.Rx.sc_psca  = sc_psca;
            e_field.Rx.sc_pext  = sc_pext;
            e_field.Rx.phi      = phi;
            e_field.Rx.phi_coil = phi_coil;
            e_field.Rx.phi_lump = phi_lump;
        end
    end
end
function[e_field] = em_efield_shield_wvie(MREDM,gpu_flag)

    func                = MREDM.functions;
    dims                = MREDM.dimensions;
    ce                  = MREDM.emc.ce;
    EP                  = MREDM.EP;
    N_wie               = dims.N_wie;
    N_shield_sie        = dims.N_shield_sie;
    N_coil              = N_shield_sie + N_wie;
    n1                  = dims.n1;
    n2                  = dims.n2;
    n3                  = dims.n3;
    ql                  = dims.ql;
    l                   = dims.l;
    mask                = dims.mask;
    idx_solution        = dims.idx_solution;
    res                 = dims.res;
    idxS                = dims.idxS;
    max_rank_wie        = dims.max_rank_N_wie;
    max_rank_shield_sie = dims.max_rank_N_shield_sie;
    Mpow                = repmat(EP.Mc_inv(idxS),[ql 1])/ce; 
    fN                  = to_GPU(MREDM.operators.N,gpu_flag);
    Z_copper_loss       = [MREDM.SIE.Z_shield_copper_loss zeros(N_shield_sie,N_wie);
                           zeros(N_wie,N_shield_sie)      MREDM.WIE.Z_copper_loss ];
    Z_lumped_loss       = [zeros(N_shield_sie,N_shield_sie) zeros(N_shield_sie,N_wie);
                           zeros(N_wie,N_shield_sie)        MREDM.WIE.Z_lumped_loss ];

    % Compute
    for i = 1:2
        if ~isempty(MREDM.fields.JcbTx) && i == 1
            Jb      = MREDM.fields.JcbTx(N_coil+1:end,:);
            Js      = MREDM.fields.JcbTx(1:N_shield_sie,:);
            Jc      = MREDM.fields.JcbTx(N_shield_sie+1:N_coil,:);
            A_ports = size(MREDM.fields.JcbTx,2);
        elseif ~isempty(MREDM.fields.JcbRx) && i == 2      
            Jb      = MREDM.fields.JcbRx(N_coil+1:end,:);
            Js      = MREDM.fields.JcbRx(1:N_shield_sie,:);
            Jc      = MREDM.fields.JcbRx(N_shield_sie+1:N_coil,:);
            A_ports = size(MREDM.fields.JcbRx,2);
        else
            continue
        end

        e_b                 = zeros(n1,n2,n3,ql,A_ports);
        pabs                = zeros(n1,n2,n3,A_ports);
        sc_pabs             = zeros(A_ports,1);
        sc_pext             = zeros(A_ports,1);
        sc_psca             = zeros(A_ports,1);


        % Total Field and Power        
        for port = 1:A_ports
            tested_einc_wire      = func.mvp_Zbw(MREDM.operators.Zbw_N,Jc(:,port),n1,n2,n3,N_wie,max_rank_wie); 
            tested_einc_shield    = func.mvp_Zbc(MREDM.operators.Zbs_N,Js(:,port),n1,n2,n3,N_sie,max_rank_shield_sie); 
            tested_einc           = tested_einc_wire+tested_einc_shield;
            Jb_port               = zeros(n1,n2,n3,ql);
            Jb_port(idx_solution) = Jb(:,port);
            tested_esca           = gather(1/ce * (func.mvp_N(to_GPU(Jb_port,gpu_flag),fN) - func.mvp_G(to_GPU(Jb_port,gpu_flag),res)));
            eb_c                  = mask.*func.mvp_invG(reshape(squeeze(tested_einc),n1,n2,n3,ql),res);
            eb_b                  = mask.*func.mvp_invG(tested_esca,res);
            e_b(:,:,:,:,port)     = eb_c + eb_b;
            pabs(:,:,:,port)      = 1/2 * EP.se.* sqrt(abs(e_b(:,:,:,1,port)).^2+abs(e_b(:,:,:,1+l,port)).^2+abs(e_b(:,:,:,1+2*l,port)).^2).^2;
            G_J                   = func.mvp_G(to_GPU(Jb_port,gpu_flag),res);
            Gx                    = gather(G_J(idx_solution));
            tested_einc           = ce*EP.Mc.*reshape(squeeze(tested_einc),n1,n2,n3,ql);
            tested_einc           = tested_einc(idx_solution);
            sc_pabs(port)         = 1/2 * real(Jb_port(idx_solution)'    * (Mpow .* Gx ));    
            sc_pext(port)         = -1/2 * real(Jb_port(idx_solution)'    * (Mpow .* tested_einc(:)));
            sc_psca(port)         = -1/2*real((Jb_port(:))' * tested_esca(:));
        end
        
        % Body Noise Covariance Matrix
        Gx = zeros(size(Jb));
        for port = 1:A_ports 
            y               = zeros(n1,n2,n3,ql);
            y(idx_solution) = Jb(:,port);
            G_J             = func.mvp_G(to_GPU(y,gpu_flag),res);
            Gx(:,port)      = gather(G_J(idx_solution));   
        end   
        phi = Jb' * (real(Mpow) .* Gx);
        phi = conj(phi);
    
        % Wire Noise Covariance Matrix
        phi_coil = zeros(A_ports,A_ports);
        phi_lump = zeros(A_ports,A_ports);
        J_loss = [Js;Jc];
        for port = 1:A_ports 
            phi_coil(port,port) = 1/2 * J_loss(:,port)' * (Z_copper_loss * J_loss(:,port));
            phi_lump(port,port) = 1/2 * J_loss(:,port)' * (Z_lumped_loss * J_loss(:,port));
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
function[e_field] = em_efield_svie_pfft_tt(MREDM,gpu_flag)

    func          = MREDM.functions;
    dims          = MREDM.dimensions;
    ce            = MREDM.emc.ce;
    EP            = MREDM.EP;
    N_sie         = dims.N_sie;
    N_shield_sie  = dims.N_shield_sie;
    N_coil        = N_shield_sie + N_sie;
    n1            = dims.n1;
    n2            = dims.n2;
    n3            = dims.n3;
    pfft_n1       = dims.pfft_n1;
    pfft_n2       = dims.pfft_n2;
    pfft_n3       = dims.pfft_n3;
    ql            = dims.ql;
    l             = dims.l;
    mask          = dims.mask;
    idx_solution  = dims.idx_solution;
    res           = dims.res;
    pfft_dims     = [pfft_n1 pfft_n2 pfft_n3 ql];
    idxS          = dims.idxS;
    Mpow          = repmat(EP.Mc_inv(idxS),[ql 1])/ce;
    fN            = to_GPU(MREDM.operators.pfft_N,gpu_flag);
    pfft_PS       = MREDM.operators.pfft_PS;
    pfft_Z_bc_N   = MREDM.operators.pfft_Z_bc_N;
    Z_copper_loss = [MREDM.SIE.Z_shield_copper_loss zeros(N_shield_sie,N_sie);
                     zeros(N_sie,N_shield_sie)      MREDM.SIE.Z_copper_loss ];
    Z_lumped_loss = [zeros(N_shield_sie,N_shield_sie) zeros(N_shield_sie,N_sie);
                     zeros(N_sie,N_shield_sie)        MREDM.SIE.Z_lumped_loss];
    
    % Compute
    for i = 1:2
        if ~isempty(MREDM.fields.JcbTx) && i == 1
            Jb      = MREDM.fields.JcbTx(N_coil+1:end,:);
            Jc      = MREDM.fields.JcbTx(N_shield_sie+1:N_coil,:);
            Js      = MREDM.fields.JcbTx(1:N_shield_sie,:);
            A_ports = size(MREDM.fields.JcbTx,2);
        elseif ~isempty(MREDM.fields.JcbRx) && i == 2        
            Jb      = MREDM.fields.JcbRx(N_coil+1:end,:);
            Jc      = MREDM.fields.JcbRx(N_shield_sie+1:N_coil,:);
            Js      = MREDM.fields.JcbRx(1:N_shield_sie,:);
            A_ports = size(MREDM.fields.JcbRx,2);
        else
            continue
        end

        e_b     = zeros(n1,n2,n3,ql,A_ports);
        pabs    = zeros(n1,n2,n3,A_ports);
        sc_pabs = zeros(A_ports,1);
        sc_pext = zeros(A_ports,1);
        sc_psca = zeros(A_ports,1);

        for port = 1:A_ports
            % Total Field and Power
            Jinc                      = reshape(pfft_PS(:,1:N_sie)*Jc(:,port),pfft_dims);
            Jsca                      = reshape(pfft_PS(:,N_sie+1:end)*Jb(:,port),pfft_dims);
            E_vox_inc                 = gather(1/ce * (func.mvp_N(to_GPU(Jinc,gpu_flag),fN) - func.mvp_G(to_GPU(Jinc,gpu_flag),res)));
            E_vox_sca                 = gather(1/ce * (func.mvp_N(to_GPU(Jsca,gpu_flag),fN) - func.mvp_G(to_GPU(Jsca,gpu_flag),res)));                
            E_vox_inc                 = pfft_Z_bc_N * Jc(:,port) +  pfft_PS(:,N_sie+1:end).' * E_vox_inc(:);
            E_vox_sca                 = pfft_PS(:,N_sie+1:end).' * E_vox_sca(:);
            tested_etot               = zeros(n1,n2,n3,ql);
            tested_esca               = zeros(n1,n2,n3,ql);
            tested_einc               = zeros(n1,n2,n3,ql);
            tested_einc_shield        = func.mvp_Zbs(MREDM.operators.tt_Zbs_N,Js(:,port));
            tested_esca(idx_solution) = E_vox_sca(:);
            tested_einc(idx_solution) = E_vox_inc(:) + tested_einc_shield(idx_solution);
            tested_etot(idx_solution) = tested_esca(idx_solution) + tested_einc(idx_solution);
            e_b(:,:,:,:,port)         = mask.*func.mvp_invG(tested_etot,res);
            pabs(:,:,:,port)          = 1/2 * EP.se .* sqrt(abs(e_b(:,:,:,1,port)).^2+abs(e_b(:,:,:,1+l,port)).^2+abs(e_b(:,:,:,1+2*l,port)).^2).^2;
            Jb_port                   = zeros(n1,n2,n3,ql);
            Jb_port(idx_solution)     = Jb(:,port);
            G_J                       = func.mvp_G(to_GPU(Jb_port,gpu_flag),res);
            Gx                        = gather(G_J(idx_solution));
            tested_einc               = ce*EP.Mc.*reshape(squeeze(tested_einc),n1,n2,n3,ql);
            tested_einc               = tested_einc(idx_solution);
            sc_pabs(port)             = 1/2 * real(Jb_port(idx_solution)'  * (Mpow .* Gx ));    
            sc_pext(port)             = 1/2 * real(Jb_port(idx_solution)' * (Mpow .* tested_einc(:)));
            sc_psca(port)             = -1/2*real((Jb_port(:))' * tested_esca(:));
        end
        
        % Body Noise Covariance Matrix
        Gx = zeros(size(Jb));
        for port = 1:A_ports 
            y               = zeros(n1,n2,n3,ql);
            y(idx_solution) = Jb(:,port);
            G_J             = func.mvp_G(to_GPU(y,gpu_flag),res);
            Gx(:,port)      = gather(G_J(idx_solution));   
        end   
        phi_2 = Jb' * (real(Mpow) .* Gx);
        phi = conj(phi_2);
    
        % Coil Noise Covariance Matrix
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
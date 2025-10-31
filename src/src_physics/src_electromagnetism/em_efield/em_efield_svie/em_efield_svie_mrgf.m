function[e_field] = em_efield_svie_mrgf(MREDM,Ue)

    se           = MREDM.EP.se;
    dims         = MREDM.dimensions;
    n1           = dims.n1;
    n2           = dims.n2;
    n3           = dims.n3;
    ql           = dims.ql;
    l            = dims.l;
    idxS         = dims.idxS;
    res          = dims.res;
    idx_solution = dims.idx_solution;

    % Wire-Coil-Shield
    if ~isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.coil) && ~isempty(MREDM.SIE.shield)
        N_wie         = dims.N_wie;
        N_sie         = dims.N_sie;
        N_shield_sie  = dims.N_shield_sie;   
        Z_copper_loss = [MREDM.SIE.Z_shield_copper_loss zeros(N_shield_sie,N_wie) zeros(N_shield_sie,N_sie);
                         zeros(N_wie,N_shield_sie)      MREDM.WIE.Z_copper_loss   zeros(N_wie,N_sie)       ;
                         zeros(N_sie,N_shield_sie)      zeros(N_sie,N_wie)        MREDM.SIE.Z_copper_loss ];
        Z_lumped_loss = [zeros(N_shield_sie,N_shield_sie) zeros(N_shield_sie,N_wie) zeros(N_shield_sie,N_sie);
                         zeros(N_wie,N_shield_sie)        MREDM.WIE.Z_lumped_loss   zeros(N_wie,N_sie)       ;
                         zeros(N_sie,N_shield_sie)        zeros(N_sie,N_wie)        MREDM.SIE.Z_lumped_loss ];
    % Wire-Coil   
    elseif ~isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.coil) && isempty(MREDM.SIE.shield)
        N_wie         = dims.N_wie;
        N_sie         = dims.N_sie;
        Z_copper_loss = [MREDM.WIE.Z_copper_loss zeros(N_wie,N_sie);
                         zeros(N_sie,N_wie)      MREDM.SIE.Z_copper_loss];
        Z_lumped_loss = [MREDM.WIE.Z_lumped_loss zeros(N_wie,N_sie);
                         zeros(N_sie,N_wie)      MREDM.SIE.Z_lumped_loss];
    % Coil-Shield
    elseif isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.coil) && ~isempty(MREDM.SIE.shield)
        N_shield_sie  = dims.N_shield_sie;
        N_sie         = dims.N_sie;
        Z_copper_loss = [MREDM.SIE.Z_shield_copper_loss zeros(N_shield_sie,N_sie);
                         zeros(N_sie,N_shield_sie)      MREDM.SIE.Z_copper_loss ];
        Z_lumped_loss = [zeros(N_shield_sie,N_shield_sie) zeros(N_shield_sie,N_sie);
                         zeros(N_sie,N_shield_sie)        MREDM.SIE.Z_lumped_loss];
    % Coil
    elseif isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.coil) && isempty(MREDM.SIE.shield)
        Z_copper_loss = MREDM.SIE.Z_copper_loss;
        Z_lumped_loss = MREDM.SIE.Z_lumped_loss;
    % Wire-Shield
    elseif ~isempty(MREDM.WIE.coil) && isempty(MREDM.SIE.coil) && ~isempty(MREDM.SIE.shield)
        N_wie         = dims.N_wie;
        N_shield_sie  = dims.N_shield_sie;
        Z_copper_loss = [MREDM.SIE.Z_shield_copper_loss zeros(N_shield_sie,N_wie);
                         zeros(N_wie,N_shield_sie)      MREDM.WIE.Z_copper_loss ];
        Z_lumped_loss = [zeros(N_shield_sie,N_shield_sie) zeros(N_shield_sie,N_wie);
                         zeros(N_wie,N_shield_sie)        MREDM.WIE.Z_lumped_loss ];
    % Wire
    elseif ~isempty(MREDM.WIE.coil) && isempty(MREDM.SIE.coil) && isempty(MREDM.SIE.shield)
        Z_copper_loss = MREDM.WIE.Z_copper_loss;
        Z_lumped_loss = MREDM.WIE.Z_lumped_loss;
    end

    % Compute
    for i = 1:2
        if ~isempty(MREDM.fields.JcbTx) && i == 1
            ac   = MREDM.fields.aTx;
            Jc   = MREDM.fields.JcbTx;
            A_ports = size(MREDM.fields.JcbTx,2);
        elseif ~isempty(MREDM.fields.JcbRx) && i == 2        
            ac   = MREDM.fields.aRx;
            Jc   = MREDM.fields.JcbRx;
            A_ports = size(MREDM.fields.JcbRx,2);
        else
            continue
        end
    
        % Total Field and Power
        e_b                 = zeros(n1*n2*n3*ql,A_ports);
        e_b(idx_solution,:) = Ue*ac;
        e_b                 = reshape(e_b,n1,n2,n3,ql,A_ports);
        pabs                = 1/2 * se.* sqrt(abs(e_b(:,:,:,1,:)).^2+abs(e_b(:,:,:,1+l,:)).^2+abs(e_b(:,:,:,1+2*l,:)).^2).^2;
    
        % Body Noise Covariance Matrix
        Gram = reshape(repelem(res^3./[1;12],[1,l-1],3),1,[]);
        phi  = zeros(A_ports,A_ports);
        for pwx = 1:ql
            Ephi = reshape(squeeze(e_b(:,:,:,pwx,:)),n1*n2*n3,A_ports);
            Ephi = Ephi(idxS,:);
            phi  = phi + Gram(pwx)*ctranspose(Ephi)*(se(idxS).*Ephi);
        end
        phi = conj(phi);
    
        % Coil Noise Covariance Matrix
        phi_coil = zeros(A_ports,A_ports);
        phi_lump = zeros(A_ports,A_ports);
        for port = 1:A_ports 
            phi_coil(port,port) = 1/2 * Jc(:,port)' * (Z_copper_loss * Jc(:,port));
            phi_lump(port,port) = 1/2 * Jc(:,port)' * (Z_lumped_loss * Jc(:,port));
        end
        phi_coil = conj(phi_coil);
        phi_lump = conj(phi_lump);
         % Save to struct and return
        if ~isempty(MREDM.fields.JcbTx) && i == 1
            e_field.Tx.e_b      = e_b;
            e_field.Tx.pabs     = pabs;
            e_field.Tx.phi      = phi;
            e_field.Tx.phi_coil = phi_coil;
            e_field.Tx.phi_lump = phi_lump;
        elseif ~isempty(MREDM.fields.JcbRx) && i == 2        
            e_field.Rx.e_b      = e_b;
            e_field.Rx.pabs     = pabs;
            e_field.Rx.phi      = phi;
            e_field.Rx.phi_coil = phi_coil;
            e_field.Rx.phi_lump = phi_lump;
        end
    end
end

function[BASIS] = em_ehfield_vie(MREDM,BASIS)

    % Preamble
    func            = MREDM.functions;
    dims            = MREDM.dimensions;
    N_ports         = size(BASIS.Ue,2);
    n1              = dims.n1;
    n2              = dims.n2;
    n3              = dims.n3;
    idxS            = dims.idxS;
    N_scat          = dims.N_scat;
    res             = dims.res;
    ql              = dims.ql;
    l               = dims.l;
    mask            = dims.mask;
    idx_solution    = dims.idx_solution;
    ce              = MREDM.emc.ce;
    mu0             = MREDM.emc.mu0;
    omega           = MREDM.emc.omega;
    k_B             = MREDM.emc.k_B;
    T0              = MREDM.emc.T0;
    M_0             = MREDM.emc.M0;
    N_conv          = 50;
    nmodes          = (unique(floor(logspace(0,log10(N_ports),N_conv))))';
    snr_num_scaling = res^3*omega*M_0/sqrt(4*k_B*T0);
    gpu_flag        = 0;

    % Compute EM fields
    fN = to_GPU(MREDM.operators.N,gpu_flag);
    fK = to_GPU(MREDM.operators.K,gpu_flag);
    for port = 1:N_ports
        Jb               = zeros(n1,n2,n3,ql);
        Jb(idx_solution) = BASIS.Jb(:,port);
        eb_b             = mask.*func.mvp_invG(gather(1/ce * (func.mvp_N(to_GPU(Jb,gpu_flag),fN) - func.mvp_G(to_GPU(Jb,gpu_flag),res))),res);
        hb_b             = mask.*func.mvp_invG(gather(func.mvp_K(to_GPU(Jb,gpu_flag),fK)),res);
        eb_b             = reshape(eb_b,n1*n2*n3*ql,1);
        hb_b             = reshape(hb_b,n1*n2*n3*ql,1);
        BASIS.Ue(:,port) = BASIS.Ue(:,port) + eb_b(idx_solution);
        BASIS.Ub(:,port) = BASIS.Ub(:,port) + hb_b(idx_solution);
    end

    % Compute Noise Covariance Matrix
    Mpow = repmat(MREDM.EP.Mc_inv(idxS),[ql 1])/ce; 
    Gx = zeros(size(BASIS.Jb));
    for port = 1:N_ports 
        y               = zeros(n1,n2,n3,ql);
        y(idx_solution) = BASIS.Jb(:,port);
        G_J             = func.mvp_G(to_GPU(y,gpu_flag),res);
        Gx(:,port)      = gather(G_J(idx_solution));   
    end   
    phi = BASIS.Jb' * (real(Mpow) .* Gx);
    BASIS.phi = conj(phi);
    BASIS = rmfield(BASIS,'Jb');

    % Compute ultimate metrics
    BASIS.cUITXE = zeros(N_scat,length(nmodes));
    BASIS.UITXE  = zeros(N_scat,1);
    BASIS.wUITXE = zeros(N_ports,N_scat);
    BASIS.cUISNR = zeros(N_scat,length(nmodes));
    BASIS.wUISNR = zeros(N_ports,N_scat);

    warning('off','MATLAB:eigs:AmbiguousSyntax');
    c = 1;
    for ii = nmodes'
        invPhi = inv(BASIS.phi(1:ii,1:ii));
        for jj = 1:N_scat
            b1m_temp           = mu0 * ( BASIS.Ub(jj,1:ii) - 1i*BASIS.Ub(jj+N_scat*l,1:ii) );
            b1p_temp           = mu0 * ( BASIS.Ub(jj,1:ii) + 1i*BASIS.Ub(jj+N_scat*l,1:ii) );
            Svox               = transpose(b1m_temp);                
            SphiS              = Svox'*(invPhi*Svox); %#ok<MINV>
            invSphiS           = inv(SphiS);
            snr_map            = 1/sqrt(invSphiS);
            BASIS.cUISNR(jj,c) = snr_num_scaling * snr_map;
            AA                 = invPhi * (b1p_temp'*b1p_temp); %#ok<MINV>
            BASIS.cUITXE(jj,c) = eigs(AA,1);
        end
        c = c+1;
    end

    UISNR_t = zeros(N_scat,1);
    UITXE_t = zeros(N_scat,1);
    invphimat = inv(BASIS.phi);
    for jj = 1:N_scat
        b1m_temp           = mu0 * ( BASIS.Ub(jj,:) - 1i*BASIS.Ub(jj+N_scat*l,:) );
        b1p_temp           = mu0 * ( BASIS.Ub(jj,:) + 1i*BASIS.Ub(jj+N_scat*l,:) );
        Svox               = transpose(b1m_temp);                
        SphiS              = Svox'*(invphimat*Svox); %#ok<MINV>
        invSphiS           = inv(SphiS);
        BASIS.wUISNR(:,jj) = transpose( invSphiS * Svox' * invphimat );%#ok<MINV>
        snr_map            = 1/sqrt(invSphiS);
        UISNR_t(jj,1)      = snr_num_scaling * snr_map;
        AA                 = conj(invphimat) * (b1p_temp'*b1p_temp); 
        [BASIS.wUITXE(:,jj), UITXE_t(jj,:)] = eigs(AA,1);
    end
    warning('on','MATLAB:eigs:AmbiguousSyntax');
    
    BASIS.UISNR       = zeros(n1,n2,n3);
    BASIS.UISNR(idxS) = UISNR_t;
    BASIS.UITXE       = zeros(n1,n2,n3);
    BASIS.UITXE(idxS) = UITXE_t;
    BASIS.nmodes      = nmodes;

end
function[h_field] = em_hfield_shield_svie(MREDM,gpu_flag)
    
    func                = MREDM.functions;
    dims                = MREDM.dimensions;
    mu0                 = MREDM.emc.mu0;
    n1                  = dims.n1;
    n2                  = dims.n2;
    n3                  = dims.n3;
    N_sie               = dims.N_sie;
    N_shield_sie        = dims.N_shield_sie;
    N_coil              = N_shield_sie + N_sie;
    ql                  = dims.ql;
    l                   = dims.l;
    mask                = dims.mask;
    res                 = dims.res;
    idx_solution        = dims.idx_solution;
    max_rank_sie        = dims.max_rank_K_sie;
    max_rank_shield_sie = dims.max_rank_K_shield_sie;
    fK                  = to_GPU(MREDM.operators.K,gpu_flag);
    
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
    
        h_b                 = zeros(n1,n2,n3,ql,A_ports);
        b1p                 = zeros(n1,n2,n3,A_ports);
        b1m                 = zeros(n1,n2,n3,A_ports);
        
        % Total Field
        for port = 1:A_ports
            tested_hinc_coil      = func.mvp_Zbc(MREDM.operators.Zbc_K,Jc(:,port),n1,n2,n3,N_sie,max_rank_sie); 
            tested_hinc_shield    = func.mvp_Zbc(MREDM.operators.Zbs_K,Js(:,port),n1,n2,n3,N_sie,max_rank_shield_sie); 
            tested_hinc           = tested_hinc_coil+tested_hinc_shield;
            Jp_port               = zeros(n1,n2,n3,ql);
            Jp_port(idx_solution) = Jb(:,port);
            tested_hsca           = gather(func.mvp_K(to_GPU(Jp_port,gpu_flag),fK)); 
            hb_c                  = mask.*func.mvp_invG(reshape(squeeze(tested_hinc),n1,n2,n3,ql),res);
            hb_b                  = mask.*func.mvp_invG(tested_hsca,res);
            h_b(:,:,:,:,port)     = hb_c + hb_b;
            b1p(:,:,:,port)       = mu0*(squeeze(h_b(:,:,:,1,port)) + 1i*squeeze(h_b(:,:,:,1+l,port)));
            b1m(:,:,:,port)       = mu0*(squeeze(h_b(:,:,:,1,port)) - 1i*squeeze(h_b(:,:,:,1+l,port)));
        end
    
        % Save to struct and return
        if ~isempty(MREDM.fields.JcbTx) && i == 1
            h_field.Tx.h_b  = h_b;
            h_field.Tx.b1p  = b1p;
            h_field.Tx.b1m  = b1m;
        elseif ~isempty(MREDM.fields.JcbRx) && i == 2        
            h_field.Rx.h_b  = h_b;
            h_field.Rx.b1p  = b1p;
            h_field.Rx.b1m  = b1m;
        end
    end
end
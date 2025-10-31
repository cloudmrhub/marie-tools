function[h_field] = em_hfield_wsvie_pfft(MREDM,gpu_flag)
   
    func         = MREDM.functions;
    dims         = MREDM.dimensions;
    mu0          = MREDM.emc.mu0;
    n1           = dims.n1;
    n2           = dims.n2;
    n3           = dims.n3;
    pfft_n1      = dims.pfft_n1;
    pfft_n2      = dims.pfft_n2;
    pfft_n3      = dims.pfft_n3;
    N_sie        = dims.N_sie;
    N_wie        = dims.N_wie;
    N_coil       = N_wie+N_sie;
    ql           = dims.ql;
    l            = dims.l;
    mask         = dims.mask;
    res          = dims.res;
    idx_solution = dims.idx_solution;
    pfft_dims    = [pfft_n1 pfft_n2 pfft_n3 ql];
    fK           = to_GPU(MREDM.operators.pfft_K,gpu_flag);
    pfft_PS      = MREDM.operators.pfft_PS;
    pfft_Z_bc_K  = MREDM.operators.pfft_Z_bc_K;

    % Compute
    for i = 1:2
        if ~isempty(MREDM.fields.JcbTx) && i == 1
            Jc      = MREDM.fields.JcbTx(1:N_coil,:);
            Jb      = MREDM.fields.JcbTx(N_coil+1:end,:);
            A_ports = size(MREDM.fields.JcbTx,2);
        elseif ~isempty(MREDM.fields.JcbRx) && i == 2        
            Jc      = MREDM.fields.JcbRx(1:N_coil,:);
            Jb      = MREDM.fields.JcbRx(N_coil+1:end,:);
            A_ports = size(MREDM.fields.JcbRx,2);
        else
            continue
        end

        h_b          = zeros(n1,n2,n3,ql,A_ports);
        b1p          = zeros(n1,n2,n3,A_ports);
        b1m          = zeros(n1,n2,n3,A_ports);
    
        for port = 1:A_ports
            Jinc                      = reshape(pfft_PS(:,1:N_coil)*Jc(:,port),pfft_dims);
            Jsca                      = reshape(pfft_PS(:,N_coil+1:end)*Jb(:,port),pfft_dims);
            H_vox_inc                 = gather(func.mvp_K(to_GPU(Jinc,gpu_flag),fK));
            H_vox_sca                 = gather(func.mvp_K(to_GPU(Jsca,gpu_flag),fK));
            H_vox_inc                 = pfft_Z_bc_K * Jc(:,port) +  pfft_PS(:,N_coil+1:end).' * H_vox_inc(:);
            H_vox_sca                 = pfft_PS(:,N_coil+1:end).' * H_vox_sca(:);
            tested_htot               = zeros(n1,n2,n3,ql);
            tested_hsca               = zeros(n1,n2,n3,ql);
            tested_hinc               = zeros(n1,n2,n3,ql);
            tested_hsca(idx_solution) = H_vox_sca(:);
            tested_hinc(idx_solution) = H_vox_inc(:);
            tested_htot(idx_solution) = tested_hsca(idx_solution) + tested_hinc(idx_solution);
            h_b(:,:,:,:,port)         = mask.*func.mvp_invG(tested_htot,res);
            b1p(:,:,:,port)           = mu0*(squeeze(h_b(:,:,:,1,port)) + 1i*squeeze(h_b(:,:,:,1+l,port)));
            b1m(:,:,:,port)           = mu0*(squeeze(h_b(:,:,:,1,port)) - 1i*squeeze(h_b(:,:,:,1+l,port)));
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
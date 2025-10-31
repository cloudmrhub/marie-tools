function[fields] = em_hfield_vie_excitation(MREDM,gpu_flag)

    % Preamble
    func            = MREDM.functions;
    dims            = MREDM.dimensions;
    N_ports         = size(MREDM.fields.Jb,2);
    n1              = dims.n1;
    n2              = dims.n2;
    n3              = dims.n3;
    res             = dims.res;
    ql              = dims.ql;
    mask            = dims.mask;
    idx_solution    = dims.idx_solution;

    % Compute H fields
    fK = to_GPU(MREDM.operators.K,gpu_flag);
    for port = 1:N_ports
        Jb               = zeros(n1,n2,n3,ql);
        Jb(idx_solution) = MREDM.fields.Jb(:,port);
        hb_b             = mask.*func.mvp_invG(gather(func.mvp_K(to_GPU(Jb,gpu_flag),fK)),res);
        hb_b             = reshape(hb_b,n1*n2*n3*ql,1);
        fields.Ub(:,port) = MREDM.BASIS.Ub(:,port) + hb_b(idx_solution);
    end

end
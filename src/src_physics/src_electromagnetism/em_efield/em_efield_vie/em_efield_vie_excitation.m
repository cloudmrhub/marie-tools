function[fields] = em_efield_vie_excitation(MREDM,gpu_flag)

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
    ce              = MREDM.emc.ce;

    % Compute E fields
    fN = to_GPU(MREDM.operators.N,gpu_flag);
    for port = 1:N_ports
        Jb                = zeros(n1,n2,n3,ql);
        Jb(idx_solution)  = MREDM.fields.Jb(:,port);
        eb_b              = mask.*func.mvp_invG(gather(1/ce * (func.mvp_N(to_GPU(Jb,gpu_flag),fN) - func.mvp_G(to_GPU(Jb,gpu_flag),res))),res);
        eb_b              = reshape(eb_b,n1*n2*n3*ql,1);
        fields.Ue(:,port) = MREDM.BASIS.Ue(:,port) + eb_b(idx_solution);
    end

end
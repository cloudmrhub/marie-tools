function[Jb,U_hat_inv,S_hat_inv,V_hat_inv] = ie_solver_vie(MREDM,gpu_flag)

    dims         = MREDM.dimensions;
    func         = MREDM.functions;
    EP           = MREDM.EP;
    N_ports      = size(MREDM.BASIS.Ue,2);
    N_b          = dims.N_b;
    ql           = dims.ql;
    idxS         = dims.idxS;
    tol          = MREDM.inputs.tol;
    ce           = MREDM.emc.ce;
    Jb           = zeros(N_b,N_ports);

    fN  = to_GPU(MREDM.operators.N,1);
    mvp = @(J)mvp_vie(J, EP, func, dims, fN);
    for i = 1:N_ports
        tic
        rhs = ce * repmat((1./EP.Mr(idxS)).* EP.Mc(idxS),ql,1).* MREDM.BASIS.Ue(:,i);
        [Jb(:,i),~,resvec] = is_gmres_vie(mvp,to_GPU(rhs,gpu_flag),50,tol,200,to_GPU(Jb(:,i),gpu_flag));
        fprintf('\tPort: %i: \n\t\tIterations: %i  \n\t\tTime: %s \n',i,size(resvec,1),hh_mm_ss(toc));
    end

    if isfield(MREDM.BASIS,'X')
        Zbb_hat_inv = MREDM.BASIS.X.'*((MREDM.BASIS.Ue.'*Jb)*MREDM.BASIS.X);
        [U_hat_inv,S_hat_inv,V_hat_inv] = svd(Zbb_hat_inv,'econ');
    else
        U_hat_inv = [];
        S_hat_inv = [];
        V_hat_inv = [];
    end

end

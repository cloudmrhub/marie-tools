function[BASIS] = ie_solver_vie_basis(MREDM,BASIS)

    dims         = MREDM.dimensions;
    func         = MREDM.functions;
    EP           = MREDM.EP;
    N_ports      = size(BASIS.Ue,2);
    N_b          = dims.N_b;
    ql           = dims.ql;
    idxS         = dims.idxS;
    tol          = MREDM.inputs.tol;
    ce           = MREDM.emc.ce;
    Jb           = zeros(N_b,N_ports);

    try 
        gpu_flag = 0;
        fN  = to_GPU(MREDM.operators.N,gpu_flag);
        mvp = @(J)mvp_vie(J, EP, func, dims, fN);
        for i = 1:N_ports
            tic
            rhs = ce * repmat((1./EP.Mr(idxS)).* EP.Mc(idxS),ql,1).* BASIS.Ue(:,i);
            [Jb(:,i),~,resvec] = is_gmres_vie(mvp,to_GPU(rhs,gpu_flag),50,tol,200,to_GPU(Jb(:,i),gpu_flag));
            fprintf('\tPort: %i: \n\t\tIterations: %i  \n\t\tTime: %s \n',i,size(resvec,1),hh_mm_ss(toc));
        end
    catch
        gpu_flag = 5;
        warning('Out of GPU memory. Running in CPU.');
        fN  = to_GPU(MREDM.operators.N,gpu_flag);
        mvp = @(J)mvp_vie(J, EP, func, dims, fN);
        for i = 1:N_ports
            tic
            rhs = ce * repmat((1./EP.Mr(idxS)).* EP.Mc(idxS),ql,1).* BASIS.Ue(:,i);
            [Jb(:,i),~,resvec] = is_gmres_vie(mvp,to_GPU(rhs,gpu_flag),50,tol,200,to_GPU(Jb(:,i),gpu_flag));
            fprintf('\tPort: %i: \n\t\tIterations: %i  \n\t\tTime: %s \n',i,size(resvec,1),hh_mm_ss(toc));
        end
    end

    if isfield(BASIS,'X')
        Zbb_hat_inv = BASIS.X.'*((BASIS.Ue.'*Jb)*BASIS.X);
        [U_hat_inv,S_hat_inv,V_hat_inv] = svd(Zbb_hat_inv,'econ');
    else
        U_hat_inv = [];
        S_hat_inv = [];
        V_hat_inv = [];
    end

    BASIS.Jb = Jb;
    BASIS.U_hat_inv = U_hat_inv;
    BASIS.S_hat_inv = S_hat_inv;
    BASIS.V_hat_inv = V_hat_inv;

end

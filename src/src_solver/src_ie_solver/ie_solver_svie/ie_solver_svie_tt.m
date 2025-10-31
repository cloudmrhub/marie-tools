function[Jcb] = ie_solver_svie_tt(MREDM,gpu_flag)

    func = MREDM.functions;
    dims = MREDM.dimensions;
    tol  = MREDM.inputs.tol;
    emc  = MREDM.emc;
    rhs  = MREDM.solver.rhs;

    A_ports      = dims.N_ports_shield;
    N_b          = dims.N_b;
    N_shield_sie = dims.N_shield_sie;
    res          = dims.res;
    idx_solution = dims.idx_solution;
    ql           = dims.ql;
    Pprec        = MREDM.solver.precond_lu.P;
    Jcb          = zeros(N_shield_sie+N_b,A_ports);

    Lprec        = to_GPU(MREDM.solver.precond_lu.L,gpu_flag);
    Uprec        = to_GPU(MREDM.solver.precond_lu.U,gpu_flag);
    Mb_prec      = to_GPU(MREDM.solver.precond_VIE.M,gpu_flag);
    N            = to_GPU(MREDM.operators.N,gpu_flag);
    Zss          = to_GPU(MREDM.SIE.Z_shield,gpu_flag);
    tt_Zbs_N     = MREDM.operators.tt_Zbs_N;
    for ij = 1:ql
        tt_Zbs_N(ij).TT.core = to_GPU(tt_Zbs_N(ij).TT.core,gpu_flag);
    end
    mvp = @(J)mvp_svie_tt(J, tt_Zbs_N, Zss, to_GPU(MREDM.EP.Mcr_inv,gpu_flag), func, N, emc, N_shield_sie, res, dims.n1, dims.n2, dims.n3, ql, idx_solution);
    for i = 1:A_ports
        tic
        [Jcb(:,i),~,resvec] = is_gmres_svie(mvp, to_GPU(rhs(:,i),gpu_flag), 50, tol, 200, Lprec, Pprec, Uprec, Mb_prec, N_shield_sie, to_GPU(Jcb(:,i),gpu_flag));
        fprintf('\t\tPort: %i, Iterations: %i, Time: %s \n',i,size(resvec,1),hh_mm_ss(toc));
    end
    
end

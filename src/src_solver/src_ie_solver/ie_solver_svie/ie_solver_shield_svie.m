function[Jcb] = ie_solver_shield_svie(MREDM,gpu_flag)

    func = MREDM.functions;
    dims = MREDM.dimensions;
    tol  = MREDM.inputs.tol;
    emc  = MREDM.emc;
    rhs  = MREDM.solver.rhs;

    N_ports             = dims.N_ports;
    N_ports_shield      = dims.N_ports_shield;
    A_ports             = N_ports_shield+N_ports;
    N_b                 = dims.N_b;
    N_sie               = dims.N_sie;
    max_rank_sie        = dims.max_rank_N_sie;
    N_shield_sie        = dims.N_shield_sie;
    max_rank_shield_sie = dims.max_rank_N_shield_sie;
    res                 = dims.res;
    n1                  = dims.n1;
    n2                  = dims.n2;
    n3                  = dims.n3;
    ql                  = dims.ql;
    idx_solution        = dims.idx_solution;
    Pprec               = MREDM.solver.precond_lu.P;
    Jcb                 = zeros(N_sie+N_b,A_ports);
    Zbs                 = MREDM.operators.Zbs_N;
    Zbc                 = MREDM.operators.Zbc_N;
    
    Lprec    = to_GPU(MREDM.solver.precond_lu.L,gpu_flag);
    Uprec    = to_GPU(MREDM.solver.precond_lu.U,gpu_flag);
    Mb_prec  = to_GPU(MREDM.solver.precond_VIE.M,gpu_flag);
    Zss      = to_GPU(MREDM.SIE.Z_shield,gpu_flag);
    Zcc      = to_GPU(MREDM.SIE.Z,gpu_flag);
    ZscU     = to_GPU(MREDM.operators.Zsc.U,gpu_flag);
    ZscV     = to_GPU(MREDM.operators.Zsc.V',gpu_flag);
    N        = to_GPU(MREDM.operators.N,gpu_flag);
    mvp      = @(J)mvp_shield_svie(J, ZscU, ZscV, Zbs, Zbc, Zss, Zcc, MREDM.EP.Mcr_inv, func, N, emc, max_rank_shield_sie, max_rank_sie, N_shield_sie, N_sie, res, n1, n2, n3, ql, idx_solution);
    for i = 1:A_ports
        tic
        [Jcb(:,i),~,resvec] = is_gmres_svie(mvp, to_GPU(rhs(:,i),gpu_flag), 50, tol, 200, Lprec, Pprec, Uprec, Mb_prec, N_shield_sie+N_sie, to_GPU(Jcb(:,i),gpu_flag));
        fprintf('\t\tPort: %i, Iterations: %i, Time: %s \n',i,size(resvec,1),hh_mm_ss(toc));
    end
    
end

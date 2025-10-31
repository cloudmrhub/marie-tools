function[Jcb] = ie_solver_wvie_pfft_tt(MREDM,gpu_flag)

    func = MREDM.functions;
    dims = MREDM.dimensions;
    tol  = MREDM.inputs.tol;
    emc  = MREDM.emc;
    rhs  = MREDM.solver.rhs;

    M_ports          = dims.M_ports;
    N_ports_shield   = dims.N_ports_shield;
    A_ports          = N_ports_shield+M_ports;
    N_b              = dims.N_b;
    N_wie            = dims.N_wie;
    N_shield_sie     = dims.N_shield_sie;
    res              = dims.res;
    idx_solution     = dims.idx_solution;
    ql               = dims.ql;
    nv               = dims.n1*dims.n2*dims.n3;
    pfft_dims        = [dims.pfft_n1 dims.pfft_n2 dims.pfft_n3 dims.ql];
    Pprec            = MREDM.solver.precond_lu.P;
    Jcb              = zeros(N_wie+N_b,A_ports);

    Lprec        = to_GPU(MREDM.solver.precond_lu.L,gpu_flag);
    Uprec        = to_GPU(MREDM.solver.precond_lu.U,gpu_flag);
    Mb_prec      = to_GPU(MREDM.solver.precond_VIE.M,gpu_flag);
    pfft_Z_cc    = to_GPU(MREDM.operators.pfft_Z_cc,gpu_flag);
    pfft_Z_bc_N  = to_GPU(MREDM.operators.pfft_Z_bc_N,gpu_flag);
    pfft_PS      = to_GPU(MREDM.operators.pfft_PS,gpu_flag);
    N            = to_GPU(MREDM.operators.pfft_N,gpu_flag);
    Multiplier_1 = to_GPU([ones(N_wie,1); -ones(N_b,1)],gpu_flag);
    Multiplier_2 = to_GPU([zeros(N_wie,1); ones(N_b,1)],gpu_flag);
    ZswU         = to_GPU(MREDM.operators.Zsw.U,gpu_flag);
    ZswV         = to_GPU(MREDM.operators.Zsw.V',gpu_flag);
    Zss          = to_GPU(MREDM.SIE.Z_shield,gpu_flag);
    tt_Zbs_N     = MREDM.operators.tt_Zbs_N;
    for ij = 1:ql
        tt_Zbs_N(ij).TT.core = to_GPU(tt_Zbs_N(ij).TT.core,gpu_flag);
    end
    mvp          = @(J)mvp_svie_pfft_tt(J, pfft_PS, pfft_Z_cc, pfft_Z_bc_N, tt_Zbs_N, Zss, ZswU, ZswV, MREDM.EP.ext_Mc_inv, func, N, Multiplier_1, Multiplier_2, emc, N_shield_sie, N_wie, res, idx_solution, pfft_dims, nv, ql);
    for i = 1:A_ports
        tic
        [Jcb(:,i),~,resvec] = is_gmres_svie(mvp, to_GPU(rhs(:,i),gpu_flag), 50, tol, 200, Lprec, Pprec, Uprec, Mb_prec, N_shield_sie+N_wie, to_GPU(Jcb(:,i),gpu_flag));
        fprintf('\t\tPort: %i, Iterations: %i, Time: %s \n',i,size(resvec,1),hh_mm_ss(toc));
    end
    
end

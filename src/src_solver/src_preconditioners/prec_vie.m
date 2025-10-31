function [MREDM] = prec_vie(MREDM)

    dimensions       = MREDM.dimensions;
    ql               = dimensions.ql;
    res              = dimensions.res;
    idx_solution     = dimensions.idx_solution;
    
    Jtemp                      = repmat(MREDM.EP.Mcr_inv, 1, 1, 1, ql);
    Mb                         = 1 / MREDM.emc.ce * MREDM.functions.mvp_G(Jtemp,res);
    MREDM.solver.precond_VIE.M = reshape(1./Mb(idx_solution),[],1);
    
end
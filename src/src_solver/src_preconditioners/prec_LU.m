function [MREDM] = prec_LU(MREDM,A)
    
    [L,U,P] = lu(A,'vector');
    MREDM.solver.precond_lu.L = inv(L);
    MREDM.solver.precond_lu.U = inv(U); 
    MREDM.solver.precond_lu.P = P; 
    
end
function [MREDM] = prec_Block_Jacobi_LU(MREDM,gpu_flag)

    lims  = [0;MREDM.WIE.coil.loop_end];
    coils = numel(lims)-1;

    Lprec = repmat(struct('L',[]), coils, 1);
    Uprec = repmat(struct('U',[]), coils, 1);
    Pprec = repmat(struct('P',[]), coils, 1);

    for k = 1:coils
        Zkk = MREDM.WIE.Z(lims(k)+1:lims(k+1),lims(k)+1:lims(k+1));                
        [Lk,Uk,Pk] = lu(Zkk,'vector');      
        Ik = eye(size(Zkk), 'like', Zkk);
        Lprec(k).L = Lk \ Ik;                
        Uprec(k).U = Ik / Uk;                
        Pprec(k).P = Pk;                   
    end

    [Lprec] = to_GPU(Lprec,gpu_flag);
    [Uprec] = to_GPU(Uprec,gpu_flag);


    MREDM.solver.precond_lu.L = @(v) prec_Block_Jacobi_LU_mult(Lprec, Pprec, Uprec, lims, v); 
    MREDM.solver.precond_lu.U = [];
    MREDM.solver.precond_lu.P = [];  
    
end
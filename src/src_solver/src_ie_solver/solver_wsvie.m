function[MREDM] = solver_wsvie(MREDM)
    
    if ~(isempty(MREDM.WIE.coil) && isempty(MREDM.SIE.coil) && isempty(MREDM.SIE.shield))
        fprintf('\tGenerate excitation ...\n');
        MREDM.solver.rhs = rhs_assembly(MREDM,MREDM.dimensions.N_scat*MREDM.dimensions.ql);
        fprintf('\tAssembling WSVIE preconditioner ...\n');
        MREDM = prec_wsvie(MREDM);
    end

    try 
        gpu_flag = 0;
        MREDM.fields.Jcb = MREDM.functions.ie_solver_wsvie(MREDM,gpu_flag);
    catch
        try
            gpu_flag = 1;
            warning('Out of GPU memory. Re-Running partially on GPU.');
            MREDM.fields.Jcb = MREDM.functions.ie_solver_wsvie(MREDM,gpu_flag);
        catch
            gpu_flag = 5;
            warning('Out of GPU memory. Running in CPU.');
            MREDM.fields.Jcb = MREDM.functions.ie_solver_wsvie(MREDM,gpu_flag);
        end
    end

    fprintf('\tComputing network parameters ...\n');
    MREDM = np_compute(MREDM);
    
end

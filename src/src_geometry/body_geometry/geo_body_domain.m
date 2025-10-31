function [MREDM] = geo_body_domain(MREDM)

    inp = MREDM.inputs;
    
    temp = load(MREDM.inputs.rhbm);    
    RHBM = temp.RHBM;
    
    MREDM.dimensions.r = RHBM.r;
    MREDM.dimensions.n1 = size(RHBM.r,1);
    MREDM.dimensions.n2 = size(RHBM.r,2);
    MREDM.dimensions.n3 = size(RHBM.r,3);
    MREDM.dimensions.res = abs(MREDM.dimensions.r(2,1,1,1) - MREDM.dimensions.r(1,1,1,1));
    MREDM.dimensions.idxS = RHBM.idxS;
    MREDM.dimensions.N_scat = size(MREDM.dimensions.idxS,1);
    MREDM.dimensions.N_b = MREDM.dimensions.N_scat*MREDM.dimensions.ql;
    mask = false(MREDM.dimensions.n1,MREDM.dimensions.n2,MREDM.dimensions.n3);
    mask(MREDM.dimensions.idxS) = true; 
    MREDM.dimensions.mask = mask;
    
    nS = length(MREDM.dimensions.idxS);
    nD = MREDM.dimensions.n1*MREDM.dimensions.n2*MREDM.dimensions.n3;
    if inp.PWX == 1
        idx_solution = zeros(12*nS,1);
        for i = 1:12
            idx_solution((i-1)*nS+1:i*nS) = (i-1)*nD + MREDM.dimensions.idxS;
        end
        idx_domain = (1:12*nD)';
    else
        idx_solution = zeros(3*nS,1);
        for i = 1:3
            idx_solution((i-1)*nS+1:i*nS) = (i-1)*nD + MREDM.dimensions.idxS;
        end
        idx_domain = (1:3*nD)';
    end
    MREDM.dimensions.idx_solution = idx_solution;
    MREDM.dimensions.idx_domain = idx_domain;

end

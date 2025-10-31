function [MREDM] = vie_assembly(MREDM,geo_flag)
     
    N       = [];
    K       = [];
    pfft_N  = [];
    pfft_K  = [];
    pfft_n  = [];
    pfft_k  = [];
    basis_N = [];
    basis_K = [];
    
    % Assemble Green's function operators
    if geo_flag == 1 && MREDM.inputs.pFFT_flag 
        if isempty(MREDM.WIE.coil) && isempty(MREDM.SIE.coil)
            n1 = MREDM.dimensions.n1;
            n2 = MREDM.dimensions.n2;
            n3 = MREDM.dimensions.n3;
            N  = assembly_N(MREDM,n1,n2,n3); 
            N  = MREDM.functions.fft_circ(N,MREDM,n1,n2,n3);
            K  = assembly_K(MREDM,n1,n2,n3); 
            K  = MREDM.functions.fft_circ(K,MREDM,n1,n2,n3);
        elseif ~isempty(MREDM.WIE.coil) || ~isempty(MREDM.SIE.coil)
            pfft_n1        = MREDM.dimensions.pfft_n1;
            pfft_n2        = MREDM.dimensions.pfft_n2;
            pfft_n3        = MREDM.dimensions.pfft_n3;
            pfft_kernel_n1 = MREDM.dimensions.pfft_kernel_n1;
            pfft_kernel_n2 = MREDM.dimensions.pfft_kernel_n1;
            pfft_kernel_n3 = MREDM.dimensions.pfft_kernel_n1;
            pfft_N         = assembly_N(MREDM,pfft_n1,pfft_n2,pfft_n3);
            pfft_n         = pfft_N(1:pfft_kernel_n1,1:pfft_kernel_n2,1:pfft_kernel_n3,:,:);
            pfft_N         = MREDM.functions.fft_circ(pfft_N,MREDM,pfft_n1,pfft_n2,pfft_n3);
            pfft_n         = MREDM.functions.fft_circ(pfft_n,MREDM,pfft_kernel_n1,pfft_kernel_n2,pfft_kernel_n3);   
            pfft_K         = assembly_K(MREDM,pfft_n1,pfft_n2,pfft_n3);
            pfft_k         = pfft_K(1:pfft_kernel_n1,1:pfft_kernel_n2,1:pfft_kernel_n3,:,:);
            pfft_K         = MREDM.functions.fft_circ(pfft_K,MREDM,pfft_n1,pfft_n2,pfft_n3);
            pfft_k         = MREDM.functions.fft_circ(pfft_k,MREDM,pfft_kernel_n1,pfft_kernel_n2,pfft_kernel_n3);
        end
    elseif geo_flag == 2 || geo_flag == 3
        n1 = MREDM.dimensions.n1;
        n2 = MREDM.dimensions.n2;
        n3 = MREDM.dimensions.n3;
        N  = assembly_N(MREDM,n1,n2,n3); 
        N  = MREDM.functions.fft_circ(N,MREDM,n1,n2,n3);
        K  = assembly_K(MREDM,n1,n2,n3); 
        K  = MREDM.functions.fft_circ(K,MREDM,n1,n2,n3);
    end

    if geo_flag == 3
        nb1 = MREDM.dimensions.nb1;
        nb2 = MREDM.dimensions.nb2;
        nb3 = MREDM.dimensions.nb3;
        basis_N = assembly_N(MREDM,nb1,nb2,nb3); 
        basis_N = MREDM.functions.fft_circ(basis_N,MREDM,nb1,nb2,nb3);
        basis_K = assembly_K(MREDM,nb1,nb2,nb3); 
        basis_K = MREDM.functions.fft_circ(basis_K,MREDM,nb1,nb2,nb3);
    end

    MREDM.operators.N = N;
    MREDM.operators.K = K;
    MREDM.operators.basis_N = basis_N;
    MREDM.operators.basis_K = basis_K;
    MREDM.operators.pfft_N = pfft_N;
    MREDM.operators.pfft_K = pfft_K;
    MREDM.operators.pfft_n = pfft_n;
    MREDM.operators.pfft_k = pfft_k;
    
end


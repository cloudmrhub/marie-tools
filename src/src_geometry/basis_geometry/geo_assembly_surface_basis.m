function [MREDM] = geo_assembly_surface_basis(MREDM)

    inp = MREDM.inputs;
    
    % Dimensions
    MREDM.dimensions.l   = 3*inp.PWX+1;
    MREDM.dimensions.ql  = 3*MREDM.dimensions.l;
    MREDM.dimensions.pwx = 9*inp.PWX+1;
    
    % Domain
    MREDM = geo_body_domain(MREDM);
    
    % Surface Basis
    coil_file              = fullfile('./data/coils/basis_files',inp.basis_support);
    coil_lumped_elements   = geo_scoil_lumped_elements(coil_file,inp.tmd);
    MREDM.SIE.coil         = geo_scoil(coil_file,coil_lumped_elements);
    MREDM.SIE.RLC          = coil_lumped_elements; 
    MREDM.dimensions.N_sie = max(MREDM.SIE.coil.index);

    % Form EP tensors
    EP = em_assembly(MREDM);

end

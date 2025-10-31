function [MREDM,EP] = geo_assembly_ultimate_basis(MREDM)

    inp = MREDM.inputs;
    
    % Dimensions
    MREDM.dimensions.l         = 3*inp.PWX+1;
    MREDM.dimensions.ql        = 3*MREDM.dimensions.l;
    MREDM.dimensions.pwx       = 9*inp.PWX+1;
    
    % Domain
    MREDM = geo_body_domain(MREDM);
    
    % Dipole Basis
    MREDM = geo_ultimate_basis(MREDM);
    % MREDM = geo_spherical_basis(MREDM);

    % Form EP tensors
    EP = em_assembly(MREDM);

end

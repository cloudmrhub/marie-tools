function [BASIS] = IncidentField_Basis(MREDM,geo_flag) 

    if geo_flag == 2
        BASIS = sSVD_rwgBasis(MREDM); 
    elseif geo_flag == 3
        BASIS = rSVD_dipoleBasis(MREDM); 
    end


end
function [MREDM,N_coil] = geo_recompute_domain(MREDM,geo_flag)
    
    inp = MREDM.inputs;
    
    MREDM = geo_body_domain(MREDM);
    
    if ~isempty(inp.coil) && isempty(inp.wire)
        if inp.pFFT_flag && geo_flag == 1 
            MREDM = pfft_surface_domain(MREDM);
            N_coil = MREDM.dimensions.N_sie;
        end
    elseif ~isempty(inp.wire) && isempty(inp.coil)
        if inp.pFFT_flag && geo_flag == 1
            MREDM = pfft_wire_domain(MREDM);
            N_coil = MREDM.dimensions.N_wie;
        end
    elseif ~isempty(inp.wire) && ~isempty(inp.coil)
        if inp.pFFT_flag && geo_flag == 1   
            MREDM = pfft_wire_surface_domain(MREDM);
            N_coil = MREDM.dimensions.N_wie + MREDM.dimensions.N_sie;
        end
    end

    MREDM.EP = em_assembly(MREDM,geo_flag);
    
end
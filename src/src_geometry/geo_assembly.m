function [MREDM] = geo_assembly(MREDM,geo_flag)
    
    % geo_flag: 1 for MARIE, 2 for SSVD, 3 for RSVD, 4 for MRGF
    inp = MREDM.inputs;
    
    MREDM.dimensions.l   = 3*inp.PWX+1;
    MREDM.dimensions.ql  = 3*MREDM.dimensions.l;
    MREDM.dimensions.pwx = 9*inp.PWX+1;
    
    MREDM = geo_body_domain(MREDM);
    
    if geo_flag == 1 || geo_flag == 4 
        if ~isempty(inp.shield)
            shield_file                   = fullfile('./data/coils/shield_files',inp.shield);
            if isfile(strcat(shield_file(1:end-4),'.json'))
                shield_lumped_elements    = geo_scoil_lumped_elements(shield_file,inp.tmd);
                MREDM.SIE.shield          = geo_scoil(shield_file,shield_lumped_elements);
                MREDM.SIE.shield_RLC      = shield_lumped_elements; 
                MREDM.SIE.coil            = [];
                MREDM.SIE.RLC             = [];
                MREDM.WIE.coil            = [];
                MREDM.WIE.RLC             = [];
            else
                MREDM.SIE.shield          = geo_shield(shield_file);
                MREDM.SIE.shield_RLC      = [];
            end
            MREDM.dimensions.N_shield_sie = max(MREDM.SIE.shield.index);
        else
            MREDM.SIE.shield              = [];
            MREDM.SIE.shield_RLC          = [];
        end
        if ~isempty(inp.coil) && isempty(inp.wire)
            coil_file                     = fullfile('./data/coils/coil_files',inp.coil);
            coil_lumped_elements          = geo_scoil_lumped_elements(coil_file,inp.tmd);
            MREDM.SIE.coil                = geo_scoil(coil_file,coil_lumped_elements);
            MREDM.SIE.RLC                 = coil_lumped_elements; 
            MREDM.dimensions.N_sie        = max(MREDM.SIE.coil.index);
            MREDM.WIE.coil                = [];
            MREDM.WIE.RLC                 = [];
            if inp.pFFT_flag && geo_flag == 1 
                MREDM = pfft_surface_domain(MREDM);
            end
        elseif ~isempty(inp.wire) && isempty(inp.coil)
            wire_file                     = fullfile('./data/coils/wire_files',inp.wire);
            wire_lumped_elements          = geo_scoil_lumped_elements(wire_file,inp.tmd);
            MREDM.WIE.coil                = geo_wcoil(wire_file,wire_lumped_elements);
            MREDM.WIE.RLC                 = wire_lumped_elements; 
            MREDM.dimensions.N_wie        = size(MREDM.WIE.coil.S_point,1);
            MREDM.SIE.coil                = [];
            MREDM.SIE.RLC                 = [];
            if inp.pFFT_flag && geo_flag == 1
                MREDM = pfft_wire_domain(MREDM);
            end
        elseif ~isempty(inp.wire) && ~isempty(inp.coil)
            coil_file                     = fullfile('./data/coils/coil_files',inp.coil);
            coil_lumped_elements          = geo_scoil_lumped_elements(coil_file,inp.tmd);
            wire_file                     = fullfile('./data/coils/wire_files',inp.wire);
            wire_lumped_elements          = geo_scoil_lumped_elements(wire_file,inp.tmd);
            MREDM.SIE.coil                = geo_scoil(coil_file,coil_lumped_elements);
            MREDM.SIE.RLC                 = coil_lumped_elements; 
            MREDM.dimensions.N_sie        = max(MREDM.SIE.coil.index);
            MREDM.WIE.coil                = geo_wcoil(wire_file,wire_lumped_elements);
            MREDM.WIE.RLC                 = wire_lumped_elements; 
            MREDM.dimensions.N_wie        = size(MREDM.WIE.coil.S_point,1);
            if inp.pFFT_flag && geo_flag == 1   
                MREDM = pfft_wire_surface_domain(MREDM);
            end
        elseif isempty(inp.wire) && isempty(inp.coil) && isempty(inp.shield)
            MREDM.BASIS.Ue = h5read(MREDM.inputs.basis_file,'/BASIS/Ue');
            MREDM.BASIS.Ub = h5read(MREDM.inputs.basis_file,'/BASIS/Ub');
        end

    elseif geo_flag == 2
        basis_file             = fullfile('./data/coils/basis_files',inp.basis_support);
        MREDM.SIE.coil         = geo_shield(basis_file);
        MREDM.dimensions.N_sie = max(MREDM.SIE.coil.index);
        MREDM.SIE.RLC        = [];
        MREDM.WIE.coil       = [];
        MREDM.WIE.RLC        = [];
        MREDM.SIE.shield     = [];
        MREDM.SIE.shield_RLC = [];
    elseif geo_flag == 3
        MREDM = geo_ultimate_basis(MREDM);
        MREDM.SIE.coil       = [];
        MREDM.SIE.RLC        = [];
        MREDM.WIE.coil       = [];
        MREDM.WIE.RLC        = [];
        MREDM.SIE.shield     = [];
        MREDM.SIE.shield_RLC = [];
    end

    MREDM.EP = em_assembly(MREDM,geo_flag);
    MREDM.functions = parse_inputs(MREDM);     
    
end

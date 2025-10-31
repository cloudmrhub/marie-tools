function [MREDM] = wsvie_coupling_assembly(MREDM) 

    Zbs_N       = [];
    Zbs_K       = [];
    Zbc_N       = [];
    Zbc_K       = [];
    Zbw_N       = [];
    Zbw_K       = [];
    Zcw         = [];
    Zsw         = [];
    Zsc         = [];
    pfft_PS     = [];
    pfft_Z_bc_N = [];
    pfft_Z_bc_K = []; 
    pfft_Z_cc   = [];
    tt_Zbs_N    = [];
    tt_Zbs_K    = [];

    % Wire-Coil-Shield
    if ~isempty(MREDM.SIE.coil) && ~isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.shield)
        fprintf('\t\tAssembling WSIE-Wire/Coil Operators ...\n'); 
        Zcw = Assembly_WSIE_block_par(MREDM.WIE.coil,MREDM.SIE.coil,MREDM.emc,MREDM.inputs); 
        fprintf('\t\tAssembling WSIE-Wire/Shield Operators ...\n'); 
        Zsw = Assembly_WSIE_block_par(MREDM.WIE.coil,MREDM.SIE.shield,MREDM.emc,MREDM.inputs); 
        Zsw.U = conj(V);
        Zsw.V = conj(U);
        fprintf('\t\tAssembling SSIE-Coil/Shield Operators ...\n'); 
        Zsc = Assembly_SIE_block_par(MREDM.SIE.coil,MREDM.SIE.shield,MREDM.emc,MREDM.inputs);
    % Wire-Coil
    elseif ~isempty(MREDM.SIE.coil) && ~isempty(MREDM.WIE.coil) &&  isempty(MREDM.SIE.shield)
        fprintf('\t\tAssembling WSIE-Wire/Coil Operators ...\n');
        Zcw = Assembly_WSIE_block_par(MREDM.WIE.coil,MREDM.SIE.coil,MREDM.emc,MREDM.inputs); 
    % Coil-Shield
    elseif ~isempty(MREDM.SIE.coil) &&  isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.shield)
        fprintf('\t\tAssembling SSIE-Coil/Shield Operators ...\n'); 
        Zsc = Assembly_SIE_block_par(MREDM.SIE.coil,MREDM.SIE.shield,MREDM.emc,MREDM.inputs);
    % Wire-Shield
    elseif isempty(MREDM.SIE.coil) && ~isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.shield)
        fprintf('\t\tAssembling WSIE-Wire/Shield Operators ...\n'); 
        Zsw = Assembly_WSIE_block_par(MREDM.WIE.coil,MREDM.SIE.shield,MREDM.emc,MREDM.inputs); 
        Zsw.U = conj(V);
        Zsw.V = conj(U);
    end

    % PFFT and TT Assembly Cases
    if MREDM.inputs.pFFT_flag
        % Wire-Coil
        if ~isempty(MREDM.SIE.coil) && ~isempty(MREDM.WIE.coil)
            fprintf('\t\tAssembling WSVIE-Wire+Coil/Body Operators ...\n'); 
            Z_cc = [MREDM.WIE.Z      Zcw.U*Zcw.V';
                    (Zcw.U*Zcw.V').' MREDM.SIE.Z ];
            N_coil = MREDM.dimensions.N_wie + MREDM.dimensions.N_sie;
            [pfft_PS,MREDM] = pfft_projection_wire_surface_assembly(MREDM);
            [Z_direct_bc_N,Z_direct_bc_K] = pfft_wire_surface_assemble_direct_bc(MREDM);
            [Z_voxel_bc_N ,Z_voxel_bc_K]  = pfft_assemble_voxel_bc(MREDM,pfft_PS,N_coil);
            pfft_Z_cc                     = pfft_assemble_voxel_cc(MREDM,pfft_PS,Z_cc,N_coil);
        % Coil
        elseif ~isempty(MREDM.SIE.coil) && isempty(MREDM.WIE.coil)
            fprintf('\t\tAssembling SVIE-Coil/Body Operators ...\n'); 
            Z_cc = MREDM.SIE.Z;
            N_coil = MREDM.dimensions.N_sie;
            [pfft_PS,MREDM] = pfft_projection_surface_assembly(MREDM);
            [Z_direct_bc_N,Z_direct_bc_K] = pfft_surface_assemble_direct_bc(MREDM);
            [Z_voxel_bc_N ,Z_voxel_bc_K]  = pfft_assemble_voxel_bc(MREDM,pfft_PS,N_coil);
            pfft_Z_cc                     = pfft_assemble_voxel_cc(MREDM,pfft_PS,Z_cc,N_coil);
        % Wire
        elseif isempty(MREDM.SIE.coil) && ~isempty(MREDM.WIE.coil)
            fprintf('\t\tAssembling WVIE-Wire/Body Operators ...\n'); 
            Z_cc = MREDM.WIE.Z;
            N_coil = MREDM.dimensions.N_wie;
            [pfft_PS,MREDM] = pfft_projection_wire_assembly(MREDM);
            [Z_direct_bc_N,Z_direct_bc_K] = pfft_wire_assemble_direct_bc(MREDM);
            [Z_voxel_bc_N ,Z_voxel_bc_K]  = pfft_assemble_voxel_bc(MREDM,pfft_PS,N_coil);
            pfft_Z_cc                     = pfft_assemble_voxel_cc(MREDM,pfft_PS,Z_cc,N_coil);
        end
        if ~isempty(MREDM.SIE.coil) || ~isempty(MREDM.WIE.coil)
            pfft_Z_bc_N = sparse(Z_direct_bc_N - Z_voxel_bc_N);
            pfft_Z_bc_K = sparse(Z_direct_bc_K - Z_voxel_bc_K);
        end
        
        % Shield
        if ~isempty(MREDM.SIE.shield)
            fprintf('\t\tAssembling SVIE-Shield/Body Operators ...\n'); 
            [tt_Zbs_N,tt_Zbs_K] = tt_shield_coupling_assembly(MREDM);
        end
        
    % Tucker Assembly Cases
    else
        % Coil
        if ~isempty(MREDM.SIE.coil)
            fprintf('\t\tAssembling SVIE-Coil/Body Operators ...\n'); 
            [Zbc_N,Zbc_K,rank_min_N_sie,rank_min_K_sie] = svie_coupling_assembly_multicolumn_Tucker(MREDM,MREDM.SIE.coil,MREDM.dimensions.N_sie);
            MREDM.dimensions.max_rank_N_sie = rank_min_N_sie;
            MREDM.dimensions.max_rank_K_sie = rank_min_K_sie;
        end
        % Wire
        if ~isempty(MREDM.WIE.coil)
            fprintf('\t\tAssembling WVIE-Wire/Body Operators ...\n'); 
            [Zbw_N,Zbw_K,rank_min_N_wie,rank_min_K_wie] = wvie_coupling_assembly_multicolumn_Tucker(MREDM);
            MREDM.dimensions.max_rank_N_wie = rank_min_N_wie;
            MREDM.dimensions.max_rank_K_wie = rank_min_K_wie;
        end
        % Shield
        if ~isempty(MREDM.SIE.shield)
            fprintf('\t\tAssembling SVIE-Shield/Body Operators ...\n'); 
            [Zbs_N,Zbs_K,rank_min_N_shield_sie,rank_min_K_shield_sie] = svie_coupling_assembly_multicolumn_Tucker(MREDM,MREDM.SIE.shield,MREDM.dimensions.N_shield_sie);
            MREDM.dimensions.max_rank_N_shield_sie = rank_min_N_shield_sie;
            MREDM.dimensions.max_rank_K_shield_sie = rank_min_K_shield_sie;
        end
    end

    MREDM.operators.Zbs_N       = Zbs_N;
    MREDM.operators.Zbs_K       = Zbs_K;
    MREDM.operators.Zbc_N       = Zbc_N;
    MREDM.operators.Zbc_K       = Zbc_K;
    MREDM.operators.Zbw_N       = Zbw_N;
    MREDM.operators.Zbw_K       = Zbw_K;
    MREDM.operators.Zcw         = Zcw;
    MREDM.operators.Zsw         = Zsw;
    MREDM.operators.Zsc         = Zsc;
    MREDM.operators.pfft_PS     = pfft_PS;
    MREDM.operators.pfft_Z_bc_N = pfft_Z_bc_N;
    MREDM.operators.pfft_Z_bc_K = pfft_Z_bc_K;
    MREDM.operators.pfft_Z_cc   = pfft_Z_cc;
    MREDM.operators.tt_Zbs_N    = tt_Zbs_N;
    MREDM.operators.tt_Zbs_K    = tt_Zbs_K;

end
function [function_handles] = parse_inputs(MREDM)

    inp = MREDM.inputs;

    % Set up Functions
    if inp.PWX == 1
        % VIE Matrix assembly
        function_handles.fft_circ                   = @assembly_fft_circ_tucker_pwl;
        % Coupling matrix-vector product functions
        function_handles.mvp_Zbs                    = @mvp_coupling_cross_TT_pwl;
        function_handles.mvp_transp_Zbs             = @mvp_coupling_cross_TT_transpose_pwl;
        function_handles.mvp_Zbc                    = @mvp_coupling_Tucker_pwl;
        function_handles.mvp_transp_Zbc             = @mvp_coupling_Tucker_transpose_pwl;
        function_handles.mvp_Zbw                    = @mvp_coupling_Tucker_pwl;
        function_handles.mvp_transp_Zbw             = @mvp_coupling_Tucker_transpose_pwl;
        % VIE Matrix matrix-vector product functions
        function_handles.mvp_N                      = @mvp_N_pwl_tucker;
        function_handles.mvp_K                      = @mvp_K_pwl_tucker;
        function_handles.mvp_G                      = @mvp_G_pwl;
        function_handles.mvp_invG                   = @mvp_invG_pwl;
        function_handles.mvp_herm_adj_K             = @mvp_K_herm_adj_pwl_tucker;
        function_handles.mvp_herm_adj_N             = @mvp_N_herm_adj_pwl_tucker;
    else
        % VIE Matrix assembly
        function_handles.fft_circ                   = @assembly_fft_circ_tucker_pwc;
        % Coupling matrix-vector product functions
        function_handles.mvp_Zbs                    = @mvp_coupling_cross_TT_pwc;
        function_handles.mvp_transp_Zbs             = @mvp_coupling_cross_TT_transpose_pwc;
        function_handles.mvp_Zbc                    = @mvp_coupling_Tucker_pwc;
        function_handles.mvp_transp_Zbc             = @mvp_coupling_Tucker_transpose_pwc;
        function_handles.mvp_Zbw                    = @mvp_coupling_Tucker_pwc;
        function_handles.mvp_transp_Zbw             = @mvp_coupling_Tucker_transpose_pwc;
        % VIE Matrix matrix-vector product functions
        function_handles.mvp_N                      = @mvp_N_pwc_tucker;
        function_handles.mvp_K                      = @mvp_K_pwc_tucker;
        function_handles.mvp_G                      = @mvp_G_pwc;
        function_handles.mvp_invG                   = @mvp_invG_pwc;
        function_handles.mvp_herm_adj_K             = @mvp_K_herm_adj_pwc_tucker;
        function_handles.mvp_herm_adj_N             = @mvp_N_herm_adj_pwc_tucker;
    end

    % Matrix-Vector Product and EM field function handles
    % Wire-Coil-Shield
    if ~isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.coil) && ~isempty(MREDM.SIE.shield)
        if inp.pFFT_flag
            function_handles.ie_solver_wsvie = @ie_solver_wsvie_pfft_tt;
            function_handles.e_field_wsvie   = @em_efield_wsvie_pfft_tt;
            function_handles.h_field_wsvie   = @em_hfield_wsvie_pfft_tt;
        else
            function_handles.ie_solver_wsvie = @ie_solver_shield_wsvie;
            function_handles.e_field_wsvie   = @em_efield_shield_wsvie;
            function_handles.h_field_wsvie   = @em_hfield_shield_wsvie;
        end
    % Wire-Coil
    elseif ~isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.coil) && isempty(MREDM.SIE.shield)
        if inp.pFFT_flag
            function_handles.ie_solver_wsvie = @ie_solver_wsvie_pfft;
            function_handles.e_field_wsvie   = @em_efield_wsvie_pfft;
            function_handles.h_field_wsvie   = @em_hfield_wsvie_pfft;
        else
            function_handles.ie_solver_wsvie = @ie_solver_wsvie;
            function_handles.e_field_wsvie   = @em_efield_wsvie;
            function_handles.h_field_wsvie   = @em_hfield_wsvie;
        end
    % Coil-Shield
    elseif isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.coil) && ~isempty(MREDM.SIE.shield)
        if inp.pFFT_flag
            function_handles.ie_solver_wsvie = @ie_solver_svie_pfft_tt;
            function_handles.e_field_wsvie   = @em_efield_svie_pfft_tt;
            function_handles.h_field_wsvie   = @em_hfield_svie_pfft_tt;
        else
            function_handles.ie_solver_wsvie = @ie_solver_shield_svie;
            function_handles.e_field_wsvie   = @em_efield_shield_svie;
            function_handles.h_field_wsvie   = @em_hfield_shield_svie;
        end
    % Coil
    elseif isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.coil) && isempty(MREDM.SIE.shield)
        if inp.pFFT_flag
            function_handles.ie_solver_wsvie = @ie_solver_svie_pfft;
            function_handles.e_field_wsvie   = @em_efield_svie_pfft;
            function_handles.h_field_wsvie   = @em_hfield_svie_pfft;
        else
            function_handles.ie_solver_wsvie = @ie_solver_svie;
            function_handles.e_field_wsvie   = @em_efield_svie;
            function_handles.h_field_wsvie   = @em_hfield_svie;
        end
    % Wire-Shield
    elseif ~isempty(MREDM.WIE.coil) && isempty(MREDM.SIE.coil) && ~isempty(MREDM.SIE.shield)
        if inp.pFFT_flag
            function_handles.ie_solver_wsvie = @ie_solver_wvie_pfft_tt;
            function_handles.e_field_wsvie   = @em_efield_wvie_pfft_tt;
            function_handles.h_field_wsvie   = @em_hfield_wvie_pfft_tt;
        else
            function_handles.ie_solver_wsvie = @ie_solver_shield_wvie;
            function_handles.e_field_wsvie   = @em_efield_shield_wvie;
            function_handles.h_field_wsvie   = @em_hfield_shield_wvie;
        end
    % Wire
    elseif ~isempty(MREDM.WIE.coil) && isempty(MREDM.SIE.coil) && isempty(MREDM.SIE.shield)
        if inp.pFFT_flag
            function_handles.ie_solver_wsvie = @ie_solver_wvie_pfft;
            function_handles.e_field_wsvie   = @em_efield_wvie_pfft;
            function_handles.h_field_wsvie   = @em_hfield_wvie_pfft;
        else
            function_handles.ie_solver_wsvie = @ie_solver_wvie;
            function_handles.e_field_wsvie   = @em_efield_wvie;
            function_handles.h_field_wsvie   = @em_hfield_wvie;
        end
    % Shield
    elseif isempty(MREDM.WIE.coil) && isempty(MREDM.SIE.coil) && ~isempty(MREDM.SIE.shield)
        if inp.pFFT_flag
            function_handles.ie_solver_wsvie = @ie_solver_svie_tt;
            function_handles.e_field_wsvie   = @em_efield_svie_tt;
            function_handles.h_field_wsvie   = @em_hfield_svie_tt;
        else
            function_handles.ie_solver_wsvie = @ie_solver_shield_only_svie;
            function_handles.e_field_wsvie   = @em_efield_shield_only_svie;
            function_handles.h_field_wsvie   = @em_hfield_shield_only_svie;
        end
    % Excitation
    elseif isempty(MREDM.WIE.coil) && isempty(MREDM.SIE.coil) && isempty(MREDM.SIE.shield)
        function_handles.ie_solver_wsvie = @ie_solver_vie;
        function_handles.e_field_wsvie   = @em_efield_vie_excitation;
        function_handles.h_field_wsvie   = @em_hfield_vie_excitation;
    end

end
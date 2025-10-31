function [MREDM] = wsvie_coupling_assembly_mrgf(MREDM) 

    Zbs_N = [];
    Zbs_K = [];
    Zbc_N = [];
    Zbc_K = [];
    Zbw_N = [];
    Zbw_K = [];
    Zcw   = [];
    Zsw   = [];
    Zsc   = [];

    % Wire-Coil-Shield
    if ~isempty(MREDM.SIE.coil) &&  ~isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.shield)
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

    % Coil
    if ~isempty(MREDM.SIE.coil)
        fprintf('\t\tAssembling SVIE-Coil/Body Operators ...\n'); 
        [Zbc_N,Zbc_K] = svie_mrgf_assembly(MREDM,MREDM.SIE.coil);
    end
    % Wire
    if ~isempty(MREDM.WIE.coil)
        fprintf('\t\tAssembling WVIE-Wire/Body Operators ...\n'); 
        [Zbw_N,Zbw_K] = wvie_mrgf_assembly(MREDM);
    end
    % Shield
    if ~isempty(MREDM.SIE.shield)
        fprintf('\t\tAssembling SVIE-Shield/Body Operators ...\n'); 
        [Zbc_N,Zbc_K] = svie_mrgf_assembly(MREDM,MREDM.SIE.shield);
    end

    MREDM.operators.Zbs_N = Zbs_N;
    MREDM.operators.Zbs_K = Zbs_K;
    MREDM.operators.Zbc_N = Zbc_N;
    MREDM.operators.Zbc_K = Zbc_K;
    MREDM.operators.Zbw_N = Zbw_N;
    MREDM.operators.Zbw_K = Zbw_K;
    MREDM.operators.Zcw   = Zcw;
    MREDM.operators.Zsw   = Zsw;
    MREDM.operators.Zsc   = Zsc;

end
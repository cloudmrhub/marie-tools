function [MREDM] = prec_wsvie(MREDM)

    % Wire-Coil-Shield
    if ~isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.coil) && ~isempty(MREDM.SIE.shield)
        Zss  = MREDM.SIE.Z_shield;
        Zsw  = MREDM.operators.Zsw;      
        ZscU = MREDM.operators.Zsc.U;    
        ZscV = MREDM.operators.Zsc.V';
        Zsc  = ZscU*ZscV;
        ZcwU = MREDM.operators.Zcw.U;  
        ZcwV = MREDM.operators.Zcw.V';
        Zcw  = ZcwU*ZcwV;
        Zww  = MREDM.WIE.Z;
        Zcc  = MREDM.SIE.Z;
        A    = [Zss   Zsw   Zsc;   
                Zsw.' Zww   Zcw;  
                Zsc.' Zcw.' Zcc];    
    % Wire-Coil   
    elseif ~isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.coil) && isempty(MREDM.SIE.shield)
        ZcwU = MREDM.operators.Zcw.U;  
        ZcwV = MREDM.operators.Zcw.V';
        Zcw  = ZcwU*ZcwV;
        Zww  = MREDM.WIE.Z;
        Zcc  = MREDM.SIE.Z;
        A    = [Zww   Zcw;  
                Zcw.' Zcc];   
    % Coil-Shield
    elseif isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.coil) && ~isempty(MREDM.SIE.shield)
        Zss  = MREDM.SIE.Z_shield;
        ZscU = MREDM.operators.Zsc.U;    
        ZscV = MREDM.operators.Zsc.V';
        Zsc  = ZscU*ZscV;
        Zcc  = MREDM.SIE.Z;
        A    = [Zss   Zsc;   
                Zsc.' Zcc]; 
    % Wire-Shield
    elseif ~isempty(MREDM.WIE.coil) && isempty(MREDM.SIE.coil) && ~isempty(MREDM.SIE.shield)
        Zss  = MREDM.SIE.Z_shield;
        Zsw  = MREDM.operators.Zsw;      
        Zww  = MREDM.WIE.Z;
        A    = [Zss   Zsw;   
                Zsw.' Zww];  
    % Coil
    elseif isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.coil) && isempty(MREDM.SIE.shield)
        A  = MREDM.SIE.Z;
    % Wire
    elseif ~isempty(MREDM.WIE.coil) && isempty(MREDM.SIE.coil) && isempty(MREDM.SIE.shield)
        A  = MREDM.WIE.Z;
    % Shield
    elseif isempty(MREDM.WIE.coil) && isempty(MREDM.SIE.coil) && ~isempty(MREDM.SIE.shield)
        A  = MREDM.SIE.Z_shield;
    end

    MREDM = prec_LU(MREDM,A);

    MREDM = prec_vie(MREDM);
    
end
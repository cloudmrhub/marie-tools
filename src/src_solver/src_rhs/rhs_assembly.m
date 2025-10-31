function [rhs] = rhs_assembly(MREDM,N_scat_ql)

    % Wire-Coil-Shield
    if ~isempty(MREDM.WIE.coil)
        M_ports          = MREDM.dimensions.M_ports;
        N_wie            = MREDM.dimensions.N_wie;
        Wie_F            = MREDM.WIE.F;
    else
        M_ports          = 0;
        N_wie            = 0;
        Wie_F            = [];
    end

    if ~isempty(MREDM.SIE.coil)
        N_ports          = MREDM.dimensions.N_ports;
        N_sie            = MREDM.dimensions.N_sie;
        Sie_F            = MREDM.SIE.F;
    else
        N_ports          = 0;
        N_sie            = 0;
        Sie_F            = [];
    end

    if ~isempty(MREDM.SIE.shield)
        N_ports_shield   = MREDM.dimensions.N_ports_shield;
        N_shield_sie     = MREDM.dimensions.N_shield_sie;
        Sie_F_shield     = MREDM.SIE.F_shield;
    else
        N_ports_shield   = 0;
        N_shield_sie     = 0;
        Sie_F_shield     = [];
    end

    rows           = N_shield_sie + N_wie + N_sie + N_scat_ql;
    cols           = N_ports_shield+M_ports+N_ports;
    rhs            = zeros(rows, cols);
    
    rows           = 1:N_shield_sie;
    cols           = 1:N_ports_shield;
    rhs(rows,cols) = Sie_F_shield;  

    chunk          = N_ports_shield;
    rows           = N_shield_sie+1:N_shield_sie+N_wie;
    cols           = chunk+1:chunk+M_ports;
    rhs(rows,cols) = Wie_F;  
    
    chunk          = N_ports_shield+M_ports;
    rows           = N_shield_sie+N_wie+1:N_shield_sie+N_wie+N_sie;
    cols           = chunk+M_ports+1:chunk+N_ports;
    rhs(rows,cols) = Sie_F;  

end
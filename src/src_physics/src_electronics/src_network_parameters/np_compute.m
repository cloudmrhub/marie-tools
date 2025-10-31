function [MREDM] = np_compute(MREDM)
    
    if isempty(MREDM.WIE.coil) && isempty(MREDM.SIE.coil) && isempty(MREDM.SIE.shield)
        return
    end

    % Wire-Coil-Shield
    if ~isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.coil) && ~isempty(MREDM.SIE.shield)
        M_ports        = MREDM.dimensions.M_ports;  
        T_ports        = MREDM.dimensions.T_ports;  
        N_ports        = MREDM.dimensions.N_ports;  
        T_ports_shield = MREDM.dimensions.T_ports_shield;  
        N_ports_shield = MREDM.dimensions.N_ports_shield; 
        N_wie          = MREDM.dimensions.N_wie;
        N_shield_sie   = MREDM.dimensions.N_shield_sie;
        
        ids = [];
        for i = 1:length(MREDM.WIE.coil.port)
            if strcmp(MREDM.WIE.coil.port(i).type,'port')
                s1 = MREDM.WIE.coil.port(i).t;
                s2 = MREDM.WIE.coil.port(i).t+1;
                ids = [ids s1 s2];
            end
        end
        ids = [1:T_ports_shield N_shield_sie+ids N_shield_sie+N_wie+1:N_shield_sie+N_wie+T_ports];
        
        It = MREDM.fields.Jcb(ids,1:N_ports_shield+M_ports+N_ports);
        Ft = MREDM.solver.rhs(ids,1:N_ports_shield+M_ports+N_ports);

    % Coil-Shield
    elseif isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.coil) && ~isempty(MREDM.SIE.shield)
        T_ports        = MREDM.dimensions.T_ports;  
        N_ports        = MREDM.dimensions.N_ports;  
        T_ports_shield = MREDM.dimensions.T_ports_shield;  
        N_ports_shield = MREDM.dimensions.N_ports_shield; 
        N_shield_sie   = MREDM.dimensions.N_shield_sie; 
        
        ids = [1:T_ports_shield N_shield_sie+1:N_shield_sie+T_ports];

        It = MREDM.fields.Jcb(ids,1:N_ports_shield+N_ports);
        Ft = MREDM.solver.rhs(ids,1:N_ports_shield+N_ports);

    % Wire-Shield
    elseif ~isempty(MREDM.WIE.coil) && isempty(MREDM.SIE.coil) && ~isempty(MREDM.SIE.shield)
        M_ports        = MREDM.dimensions.M_ports;  
        T_ports_shield = MREDM.dimensions.T_ports_shield;  
        N_ports_shield = MREDM.dimensions.N_ports_shield; 
        N_shield_sie   = MREDM.dimensions.N_shield_sie; 
        
        ids = [];
        for i = 1:length(MREDM.WIE.coil.port)
            if strcmp(MREDM.WIE.coil.port(i).type,'port')
                s1 = MREDM.WIE.coil.port(i).t;
                s2 = MREDM.WIE.coil.port(i).t+1;
                ids = [ids s1 s2];
            end
        end
        ids = [1:T_ports_shield N_shield_sie+1:N_shield_sie+ids];
        
        It = MREDM.fields.Jcb(ids,1:N_ports_shield+M_ports);
        Ft = MREDM.solver.rhs(ids,1:N_ports_shield+M_ports);

    % Wire-Coil  
    elseif ~isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.coil) && isempty(MREDM.SIE.shield)
        M_ports = MREDM.dimensions.M_ports;  
        T_ports = MREDM.dimensions.T_ports;  
        N_ports = MREDM.dimensions.N_ports; 
        N_wie   = MREDM.dimensions.N_wie;
        
        ids = [];
        for i = 1:length(MREDM.WIE.coil.port)
            if strcmp(MREDM.WIE.coil.port(i).type,'port')
                s1 = MREDM.WIE.coil.port(i).t;
                s2 = MREDM.WIE.coil.port(i).t+1;
                ids = [ids s1 s2];
            end
        end
        ids = [ids N_wie+1:N_wie+T_ports];
        
        It = MREDM.fields.Jcb(ids,1:M_ports+N_ports);
        Ft = MREDM.solver.rhs(ids,1:M_ports+N_ports);

    % Coil
    elseif isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.coil) && isempty(MREDM.SIE.shield)  
        T_ports = MREDM.dimensions.T_ports;  
        N_ports = MREDM.dimensions.N_ports;  
        
        It = MREDM.fields.Jcb(1:T_ports,1:N_ports);
        Ft = MREDM.solver.rhs(1:T_ports,1:N_ports);

    % Wire
    elseif ~isempty(MREDM.WIE.coil) && isempty(MREDM.SIE.coil) && isempty(MREDM.SIE.shield)
        M_ports = MREDM.dimensions.M_ports;  
        
        ids     = [];
        for i = 1:length(MREDM.WIE.coil.port)
            if strcmp(MREDM.WIE.coil.port(i).type,'port')
                s1  = MREDM.WIE.coil.port(i).t;
                s2  = MREDM.WIE.coil.port(i).t+1;
                ids = [ids s1 s2];
            end
        end
        
        It = MREDM.fields.Jcb(ids,1:M_ports);
        Ft = MREDM.solver.rhs(ids,1:M_ports);

    % Shield
    elseif isempty(MREDM.WIE.coil) && isempty(MREDM.SIE.coil) && ~isempty(MREDM.SIE.shield)  
        T_ports_shield = MREDM.dimensions.T_ports_shield;  
        N_ports_shield = MREDM.dimensions.N_ports_shield;  
        
        It = MREDM.fields.Jcb(1:T_ports_shield,1:N_ports_shield);
        Ft = MREDM.solver.rhs(1:T_ports_shield,1:N_ports_shield);
    end

    Ip = -transpose(Ft)*It;
    MREDM.fields.netp.YP_s = (Ip + Ip.')/2;
    MREDM.fields.netp.ZP_s = np_y2z(MREDM.fields.netp.YP_s);
    MREDM.fields.netp.SP_s = np_z2s(MREDM.fields.netp.ZP_s,MREDM.emc.Z0);
    
end
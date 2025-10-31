function [MREDM] = shield_sie_assembly(MREDM)

    if ~isempty(MREDM.SIE.shield)

        % Assembly Galerkin Matrix of EFIE
        Shield = MREDM.SIE.shield;
        N_shield_sie = MREDM.dimensions.N_shield_sie;
        [Z,Z_copper_loss] = Assembly_SIE_par(Shield,MREDM.emc,MREDM.inputs);
        
        if ~isempty(MREDM.SIE.shield_RLC)
            
            Nloads = length(Shield.port);
            N_ports = 0;
            for i = 1:length(Shield.port)
                if strcmp(Shield.port(i).type,'port')
                    N_ports = N_ports+1;
                end
            end
    
            F   = zeros(max(Shield.index),N_ports);   
            % ---------------------------------------------------------------------
            %             delta-gap method 
            % ---------------------------------------------------------------------
            
            Z_lumped_loss = sparse(N_shield_sie,N_shield_sie);
            c = 1;
            for i = 1:Nloads
                if strcmp(Shield.port(i).type,'element')
                    [Z,Z_lumped_loss] = assembly_le(Z,Z_lumped_loss,Shield,Shield.port(i),MREDM.emc);
                end
                if strcmp(Shield.port(i).type,'port')
                    Vp = excitation_coil(Shield,Shield.port(i));
                    F(:,c)    = Vp;
                    c = c+1;
                end
            end
        
            idx                               = find(sum(abs(F),2));
            T_ports                           = length(idx);
            MREDM.dimensions.N_ports_shield   = N_ports;
            MREDM.dimensions.T_ports_shield   = T_ports;
            
            MREDM.SIE.F_shield = F;
            MREDM.SIE.Z_shield_lumped_loss = Z_lumped_loss;


        else
            MREDM.dimensions.N_ports_shield   = 0;
            MREDM.dimensions.T_ports_shield   = 0;
            MREDM.SIE.F_shield   = [];
            MREDM.SIE.Z_shield_lumped_loss = [];
        end

        MREDM.SIE.Z_shield = -Z;
        MREDM.SIE.Z_shield_copper_loss = Z_copper_loss;
    
    end

end
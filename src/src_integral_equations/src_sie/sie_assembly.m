function [MREDM] = sie_assembly(MREDM)

    if ~isempty(MREDM.SIE.coil)
        % Assembly Galerkin Matrix of EFIE
        Scoil = MREDM.SIE.coil;
        N_sie = MREDM.dimensions.N_sie;
        Nloads = length(Scoil.port);
        
        N_ports = 0;
        for i = 1:length(Scoil.port)
            if strcmp(Scoil.port(i).type,'port')
                N_ports = N_ports+1;
            end
        end
    
        [Z,Z_copper_loss] = Assembly_SIE_par(Scoil,MREDM.emc,MREDM.inputs);
        F   = zeros(max(Scoil.index),N_ports);   
        
        % ---------------------------------------------------------------------
        %             delta-gap method 
        % ---------------------------------------------------------------------
        
        Z_lumped_loss = sparse(N_sie,N_sie);
        c = 1;
        for i = 1:Nloads
            if strcmp(Scoil.port(i).type,'element')
                [Z,Z_lumped_loss] = assembly_le(Z,Z_lumped_loss,Scoil,Scoil.port(i),MREDM.emc);
            end
            if strcmp(Scoil.port(i).type,'port')
                Vp = excitation_coil(Scoil,Scoil.port(i));
                F(:,c)    = Vp;
                c = c+1;
            end
        end
        
        idx                        = find(sum(abs(F),2));
        T_ports                    = length(idx);
        MREDM.dimensions.N_ports   = N_ports;
        MREDM.dimensions.T_ports   = T_ports;
        
        MREDM.SIE.F = F;
        MREDM.SIE.Z = -Z;
        MREDM.SIE.Z_copper_loss = Z_copper_loss;
        MREDM.SIE.Z_lumped_loss = Z_lumped_loss;

    end

end
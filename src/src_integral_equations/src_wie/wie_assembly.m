function [MREDM] = wie_assembly(MREDM)

    if ~isempty(MREDM.WIE.coil)
        % Assembly Galerkin Matrix of WIE
        Dwire   = MREDM.emc.thick_wire;
        Wcoil   = MREDM.WIE.coil;
        N_wie   = MREDM.dimensions.N_wie;
        Ports   = Wcoil.port;
        quad_od = MREDM.inputs.Quad_order_wie;
        Nloads  = length(Ports);
        
        N_ports = 0;
        for i = 1:Nloads
            if strcmp(Ports(i).type,'port')
                N_ports = N_ports+1;
            end
        end
        
        [Z_EFIE,Z_copper_loss] = Assembly_WIE_triangle_basis_vec(Wcoil,Dwire/2,quad_od,MREDM.emc);
        
        Z_lumped_loss = sparse(N_wie,N_wie);
        for i = 1:Nloads
            if strcmp(Ports(i).type,'element')
                [Z_EFIE,Z_lumped_loss] = assembly_wle(Z_EFIE,Z_lumped_loss,Ports(i),MREDM.emc);
            end
        end
        
        F = excitation_wire(Wcoil,N_wie,N_ports,Nloads);
    
        MREDM.WIE.Z              = -Z_EFIE;
        MREDM.WIE.Z_copper_loss  = Z_copper_loss;
        MREDM.WIE.Z_lumped_loss  = Z_lumped_loss;
        MREDM.WIE.F              = F;
        MREDM.dimensions.M_ports = N_ports;

    end
    
end
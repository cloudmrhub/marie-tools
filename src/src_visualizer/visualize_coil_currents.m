function visualize_coil_currents(MREDM,port)
    
    fs = 15;

    dims   = MREDM.dimensions;
    fields = MREDM.fields;

    % Wire-Coil-Shield
    if ~isempty(MREDM.WIE.coil)
        M_ports          = dims.M_ports;
        N_wie            = dims.N_wie;
    else
        M_ports          = 0;
        N_wie            = 0;
    end
    if ~isempty(MREDM.SIE.coil)
        N_ports          = dims.N_ports;
        N_sie            = dims.N_sie;
    else
        N_ports          = 0;
        N_sie            = 0;
    end
    if ~isempty(MREDM.SIE.shield)
        N_ports_shield   = dims.N_ports_shield;
        N_shield_sie     = dims.N_shield_sie;
    else
        N_ports_shield   = 0;
        N_shield_sie     = 0;
    end
    
    for i = 1:2
        if i == 1 && isfield(fields,'JcbRx')
            Jcb = fields.JcbRx;
        elseif i == 2 && isfield(fields,'JcbTx') 
            Jcb = fields.JcbTx;
        else
            continue;
        end
        rows       = 1:N_shield_sie;
        cols       = 1:N_ports_shield;
        Jcb_shield = Jcb(rows,cols);
    
        chunk      = N_ports_shield;
        rows       = N_shield_sie+1:N_shield_sie+N_wie;
        cols       = chunk+1:chunk+M_ports;
        Jcb_wire   = Jcb(rows,cols);
        
        chunk      = N_ports_shield+M_ports;
        rows       = N_shield_sie+N_wie+1:N_shield_sie+N_wie+N_sie;
        cols       = chunk+1:chunk+N_ports;
        Jcb_coil   = Jcb(rows,cols);


        v1 = []; v2 = []; v3 = []; c1 = []; c2 = []; c3 = [];
        if isfield(COIL,'shield')
            [v1,c1] = surface_to_cartesian_currents(Jcb_shield(:,port),MREDM.SIE.shield);
        end
    
        if isfield(WIRE,'coil')
            [v2,c2] = line_to_cartesian_currents(Jcb_wire(:,port),MREDM.WIE.coil);
        end
    
        if isfield(COIL,'coil')
            [v3,c3] = surface_to_cartesian_currents(Jcb_coil(:,port),MREDM.SIE.coil);
        end

        c = [c1;c2;c3];
        v = [v1;v2;v3];

        figure
        quiver3(c(:,1),c(:,2),c(:,3),v(:,1),v(:,2),v(:,3),1,'k'); 
        xlabel('x','interpreter','latex','Fontsize',fs);
        ylabel('y','interpreter','latex','Fontsize',fs);
        zlabel('z','interpreter','latex','Fontsize',fs);
        ax = gca;
        ax.FontSize = fs;
        axis equal;
        axis on;
        grid on;

        if i == 1 && isfield(fields,'JcbRx')
            title('$Re(j_{\rm Rx})$','Interpreter','latex','Fontsize',fs);
        elseif i == 2 && isfield(fields,'JcbTx') 
            title('$Re(j_{\rm Tx})$','Interpreter','latex','Fontsize',fs);
        end

    end

end
function visualize_geometry(MREDM)

    fs = 15;

    dims = MREDM.dimensions;
    WIRE = MREDM.WIE;
    COIL = MREDM.SIE;
    
    % ==========================================================
    % Body
    % ==========================================================
    xd = dims.r(:,:,:,1);
    yd = dims.r(:,:,:,2);
    zd = dims.r(:,:,:,3);
    id = dims.idxS;
    P  = [xd(id) yd(id) zd(id)];
    k  = boundary(P,1);        
    trisurf(k,P(:,1),P(:,2),P(:,3),'LineStyle','none','FaceColor',[254, 227, 212]./255,'FaceAlpha',1); hold on;
    
    % ==========================================================
    % Wire
    % ==========================================================
    if isfield(WIRE,'coil')
        xw = [WIRE.F_point(:, 1)'; WIRE.S_point(:, 1)'; NaN(1, size(WIRE.F_point, 1))];
        yw = [WIRE.F_point(:, 2)'; WIRE.S_point(:, 2)'; NaN(1, size(WIRE.F_point, 1))];
        zw = [WIRE.F_point(:, 3)'; WIRE.S_point(:, 3)'; NaN(1, size(WIRE.F_point, 1))];
        xw = xw(:);
        yw = yw(:);
        zw = zw(:);
        plot3(xw,yw,zw,'Color',[184, 115, 51]/255,'LineStyle','-','LineWidth',2); hold on;
    end
    
    % ==========================================================
    % Surface Coil
    % ==========================================================
    if isfield(COIL,'coil')
        TR1 = triangulation(COIL.coil.elem(1:3,:)',COIL.coil.node');
        trisurf(TR1,'FaceColor',[184, 115, 51]/255,'FaceAlpha',1,'LineStyle','-');  hold on;        
    end

    % ==========================================================
    % Surface Shield
    % ==========================================================
    if isfield(COIL,'shield')
        TR2 = triangulation(COIL.shield.elem(1:3,:)',COIL.shield.node');
        trisurf(TR2,'FaceColor',[184, 115, 51]/255,'FaceAlpha',0.3,'LineStyle','none');  hold on;
    end

    % ==========================================================
    % Dipole Basis Support
    % ==========================================================
    if isfield(dims,'r_basis')
        xb = dims.r_basis(:,:,:,1);
        yb = dims.r_basis(:,:,:,2);
        zb = dims.r_basis(:,:,:,3);
        ib = dims.basis_idxS;
        scatter3(xb(ib),yb(ib),zb(ib),20,'sq','MarkerEdgeColor',[255, 165, 0]/255,'MarkerFaceColor',[255, 165, 0]/255,'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);  hold on;
    end
    
    xlabel('x (m)','interpreter','latex','Fontsize',fs);
    ylabel('y (m)','interpreter','latex','Fontsize',fs);
    zlabel('z (m)','interpreter','latex','Fontsize',fs);
    
    ax = gca;
    ax.FontSize = fs;
    axis equal;
    axis on;
    grid on;
    camlight(-37,18)
    lighting gouraud;

    hold off;
    
end
    
% pos = zeros(length(COIL.port),3);
% s = zeros(length(COIL.port),1);
% for i = 1:length(COIL.port)
%     pos(i,:) = COIL.node(:,COIL.edge(1,COIL.index == COIL.port(i).t(1)));
%     s(i) = i;
% end
% s = string(s);
% text(pos(:,1),pos(:,2),pos(:,3),s,'Interpreter','none','HorizontalAlignment','center','Color','k','FontSize',fs+5,'FontWeight','bold'); hold on;
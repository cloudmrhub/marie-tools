function visualize_body_coil_wire(r,idxS,COIL,WIRE)

    figure
    fs = 15;
    
    xd = r(:,:,:,1);
    yd = r(:,:,:,2);
    zd = r(:,:,:,3);
    id = idxS;
    
    
    P = [xd(id) yd(id) zd(id)];
    
    k = boundary(P,1);
    
    trisurf(k,P(:,1),P(:,2),P(:,3),'LineStyle','none','FaceColor',[254, 227, 212]./255,'FaceAlpha',1);
    camlight(-37,18)
    lighting gouraud;
    hold on;

    xw = [WIRE.F_point(:, 1) WIRE.S_point(:, 1)]';
    yw = [WIRE.F_point(:, 2) WIRE.S_point(:, 2)]';
    zw = [WIRE.F_point(:, 3) WIRE.S_point(:, 3)]';
    
    TR = triangulation(COIL.elem(1:3,:)',COIL.node');

    plot3(xw,yw,zw,'Color',[184, 115, 51]/255,'LineStyle','-','LineWidth',2); 
    hold on;
    
    trisurf(TR,'FaceColor',[184, 115, 51]/255,'FaceAlpha',1,'LineStyle','none'); 
    xlabel('x (m)','interpreter','latex','Fontsize',fs);
    ylabel('y (m)','interpreter','latex','Fontsize',fs);
    zlabel('z (m)','interpreter','latex','Fontsize',fs);
    ax = gca;
    ax.FontSize = fs;
    axis equal;
    axis off;
    
end

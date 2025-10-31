function visualize_body_wire(r,idxS,COIL)
    
    figure
    fs = 15;
    
    xd = r(:,:,:,1);
    yd = r(:,:,:,2);
    zd = r(:,:,:,3);
    id = idxS;
    

    P = [xd(id) yd(id) zd(id)];
    
    k = boundary(P,1);
    
    trisurf(k,P(:,1),P(:,2),P(:,3),'LineStyle','none','FaceColor',[254, 227, 212]./255,'FaceAlpha',1);
    camlight(-90,18)
    lighting gouraud;
    hold on;

    xw = [COIL.F_point(:, 1)'; COIL.S_point(:, 1)'; NaN(1, size(COIL.F_point, 1))];
    yw = [COIL.F_point(:, 2)'; COIL.S_point(:, 2)'; NaN(1, size(COIL.F_point, 1))];
    zw = [COIL.F_point(:, 3)'; COIL.S_point(:, 3)'; NaN(1, size(COIL.F_point, 1))];
    
    xw = xw(:);
    yw = yw(:);
    zw = zw(:);
    
    plot3(xw, yw, zw, 'Color', [184, 115, 51]/255, 'LineStyle', '-', 'LineWidth', 2);

    xlabel('x (m)','interpreter','latex','Fontsize',fs);
    ylabel('y (m)','interpreter','latex','Fontsize',fs);
    zlabel('z (m)','interpreter','latex','Fontsize',fs);
    ax = gca;
    ax.FontSize = fs;
    axis equal;
    axis off;
    
end
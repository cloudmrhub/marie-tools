function visualize_body_coil(r,idxS,COIL)

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
    
    % xb = r(:,:,:,1);
    % yb = r(:,:,:,2);
    % zb = r(:,:,:,3);
    % ib = idxS;
    % scatter3(xb(ib),yb(ib),zb(ib),20,'sq','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);  hold on;
    
    TR = triangulation(COIL.elem(1:3,:)',COIL.node');
    
    trisurf(TR,'FaceColor',[184, 115, 51]/255,'FaceAlpha',0.4,'LineStyle','none'); 
    xlabel('x (m)','interpreter','latex','Fontsize',fs);
    ylabel('y (m)','interpreter','latex','Fontsize',fs);
    zlabel('z (m)','interpreter','latex','Fontsize',fs);
    ax = gca;
    ax.FontSize = fs;
    axis equal;
    axis off;
    
end

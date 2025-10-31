function visualize_body_coil_shield(r,idxS,COIL,SHIELD)

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
    
    TR1 = triangulation(COIL.elem(1:3,:)',COIL.node');
    TR2 = triangulation(SHIELD.elem(1:3,:)',SHIELD.node');
    
    trisurf(TR1,'FaceColor',[184, 115, 51]/255,'FaceAlpha',1,'LineStyle','none');  hold on;
    trisurf(TR2,'FaceColor',[184, 115, 51]/255,'FaceAlpha',0.3,'LineStyle','none');  hold on;
    xlabel('x (m)','interpreter','latex','Fontsize',fs);
    ylabel('y (m)','interpreter','latex','Fontsize',fs);
    zlabel('z (m)','interpreter','latex','Fontsize',fs);
    ax = gca;
    ax.FontSize = fs;
    axis equal;
    axis off;
    
end

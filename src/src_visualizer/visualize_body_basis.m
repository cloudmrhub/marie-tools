function visualize_body_basis(MREDM)
    
    figure
    fs = 15;
    
    xd = MREDM.dimensions.r(:,:,:,1);
    yd = MREDM.dimensions.r(:,:,:,2);
    zd = MREDM.dimensions.r(:,:,:,3);
    id = MREDM.dimensions.idxS;
    Pd = [xd(id) yd(id) zd(id)];
    kd = boundary(Pd,1);
    
    xb = MREDM.dimensions.r_basis(:,:,:,1);
    yb = MREDM.dimensions.r_basis(:,:,:,2);
    zb = MREDM.dimensions.r_basis(:,:,:,3);
    ib = MREDM.dimensions.basis_idxS;

    trisurf(kd,Pd(:,1),Pd(:,2),Pd(:,3),'LineStyle','none','FaceColor',[254, 227, 212]./255,'FaceAlpha',1); hold on;
    camlight(-37,18)
    lighting gouraud;
    scatter3(xb(ib),yb(ib),zb(ib),20,'sq','MarkerEdgeColor',[255, 165, 0]/255,'MarkerFaceColor',[255, 165, 0]/255,'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);  hold on;
    xlabel('x (m)','interpreter','latex','Fontsize',fs);
    ylabel('y (m)','interpreter','latex','Fontsize',fs);
    zlabel('z (m)','interpreter','latex','Fontsize',fs);
    ax = gca;
    ax.FontSize = fs;
    axis equal;
    axis off;
    view([45 45]);
    
    % scatter3(xd(id),yd(id),zd(id),20,'sq','MarkerEdgeColor','b','MarkerFaceColor','b'); 

end


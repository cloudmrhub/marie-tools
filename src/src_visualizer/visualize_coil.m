function visualize_coil(COIL)
    
    figure
    fs = 15;
    
    TR = triangulation(COIL.elem(1:3,:)',COIL.node');
    
    trisurf(TR,'LineStyle','none','FaceColor',[0.72 0.45 0.2],'FaceAlpha',1); 
    xlabel('x (m)','interpreter','latex','Fontsize',fs);
    ylabel('y (m)','interpreter','latex','Fontsize',fs);
    zlabel('z (m)','interpreter','latex','Fontsize',fs);
    ax = gca;
    ax.FontSize = fs;
    axis equal;
    axis off;
    
end
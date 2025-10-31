function visualize_body(r,idxS)

    figure
    fs = 15;
    
    xd = r(:,:,:,1);
    yd = r(:,:,:,2);
    zd = r(:,:,:,3);
    id = idxS;
    
    
    P = [xd(id) yd(id) zd(id)];
    
    k = boundary(P,1);
    
    trisurf(k,P(:,1),P(:,2),P(:,3),'LineStyle','none','FaceColor',[254, 227, 212]./255);
    camlight
    lighting gouraud;
    hold on;
    ax = gca;
    ax.FontSize = fs;
    axis equal;
    axis off;
    
end
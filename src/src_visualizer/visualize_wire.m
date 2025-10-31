function visualize_wire(COIL)
    
    % figure
    
    fs = 15;
    
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
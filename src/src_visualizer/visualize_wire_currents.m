function visualize_wire_currents(J,wire)

    figure
    
    [v,c] = line_to_cartesian_currents(J,wire);

    quiver3(c(:,1),c(:,2),c(:,3),v(:,1),v(:,2),v(:,3),1,'k'); hold on;
    xlabel('x','interpreter','latex');
    ylabel('y','interpreter','latex');
    zlabel('z','interpreter','latex');
    axis equal;
    grid on;
    title('$Re(j_{\rm wire})$','Interpreter','latex');

end
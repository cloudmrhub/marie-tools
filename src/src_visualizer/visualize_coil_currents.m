function visualize_coil_currents(J,coil,pipi)

    % figure

    [c, j_cart] = surface_to_cartesian_currents(J,coil);

    expt = exp(1i*pipi);
    jx = imag(j_cart(:,1)*expt);
    jy = imag(j_cart(:,2)*expt);
    jz = imag(j_cart(:,3)*expt);
    v = [jx jy jz];

    quiver3(c(:,1),c(:,2),c(:,3),v(:,1),v(:,2),v(:,3),3,'k'); hold on;
    xlabel('x','interpreter','latex');
    ylabel('y','interpreter','latex');
    zlabel('z','interpreter','latex');
    axis equal;
    grid on;
    axis off;

end
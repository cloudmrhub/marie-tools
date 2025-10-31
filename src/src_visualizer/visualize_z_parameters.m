function visualize_z_parameters(ZPm)
    
    figure

    colorMap = [linspace(0,1,256)', zeros(256,2)];
    MrZ = max(max(abs(real(ZPm))));
    MiZ = max(max(abs(imag(ZPm))));
    Nports = size(ZPm,1);
    
    subplot(1,2,1)
    imagesc(abs(real(ZPm)),[0 70]);
    for i = 1:Nports
        for j = 1:Nports
            text(i-0.5,j,num2str(abs(real(ZPm(i,j))), '%2.1f'), 'Color', [1-abs(real(ZPm(i,j)))/MrZ 1 1]);
        end
    end
    colormap(colorMap);
    cor = colorbar;
    cor.Ticks = [0,10,20,30,40,50,60,70];
    set(gca,'xtick',1:Nports);
    set(gca,'ytick',1:Nports);
    xlabel('Port','interpreter','latex');
    ylabel('Port','interpreter','latex');
    title('Real$\{Z_m\}$ parameters','interpreter','latex');
    axis image;
    
    subplot(1,2,2)
    imagesc(abs(imag(ZPm)),[0 35]);
    for i = 1:Nports
        for j = 1:Nports
            text(i-0.5,j,num2str(abs(imag(ZPm(i,j))), '%2.1f'), 'Color', [1-abs(imag(ZPm(i,j)))/MiZ 1 1]);
        end
    end
    
    colormap(colorMap);
    cor = colorbar;
    cor.Ticks = [0,5,10,15,20,25,30,35];
    set(gca,'xtick',1:Nports);
    set(gca,'ytick',1:Nports);
    xlabel('Port','interpreter','latex');
    ylabel('Port','interpreter','latex');
    title('Imag$\{Z_m\}$ parameters','interpreter','latex');
    axis image;

end
function visualize_s_parameters(SPm)
    
    fs = 24;
    set(groot,'defaultAxesTickLabelInterpreter','latex');  

    colorMap = [linspace(0,1,256)', zeros(256,2)];
    MS = min(min(20*log10(abs(SPm))));
    Nports = size(SPm,1);
    
    figure
    imagesc(20*log10(abs(SPm)),[-40 -5]);
    for i = 1:Nports
        for j = 1:Nports
            splot = 20*log10(abs(SPm(i,j)));
            if isinf(splot)
                fracp = 0;
            else
                fracp = splot/MS;
            end
            text(i-0.5,j,num2str_custom(splot), 'Color', [fracp 1 1], 'FontSize', fs);
        end
    end
    colormap(colorMap);
    cor = colorbar;
    cor.Ticks = [-35,-30,-25,-20,-15,-10,-5];
    cor.Title.String = " ";
    cor.Title.FontSize = fs;
    cor.FontSize = fs;
    cor.Title.Interpreter = "latex";
    cor.TickLabelInterpreter = "latex";
    cor.Label.Interpreter = "latex";
    ax = gca;
    ax.FontSize = fs;
    set(ax,'xtick',1:Nports);
    set(ax,'ytick',1:Nports);
    xlabel('Port','interpreter','latex','Fontsize',fs);
    ylabel('Port','interpreter','latex','Fontsize',fs);
    title('$20\log_{10}(|S_m|)$ parameters','interpreter','latex','Fontsize',fs);
    axis image;

end

function[b] = num2str_custom(a)

    if isinf(a)
        b = "0";
    else
        b = num2str(a, '%2.1f');
    end

end
function visualize_ZPm_SPm_freq_sweep(figure_freq)

    ZPm_plot_Rx = figure_freq.ZPm_plot_Rx;
    SPm_plot_Rx = figure_freq.SPm_plot_Rx;
    ZPm_plot_Tx = figure_freq.ZPm_plot_Tx;
    SPm_plot_Tx = figure_freq.SPm_plot_Tx;
    freq_range  = figure_freq.freq_range;
    freq_idx    = figure_freq.freq_idx;

    for i = 1:2
        if (i == 1 && ~isempty(ZPm_plot_Tx))
            ZPm_plot = ZPm_plot_Tx;
            SPm_plot = SPm_plot_Tx;
            mode     = 'Tx';
        elseif (i == 2 && ~isempty(ZPm_plot_Rx))
            ZPm_plot = ZPm_plot_Rx;
            SPm_plot = SPm_plot_Rx;
            mode     = 'Rx';
        else
            continue
        end

        figure(i)
        hold on;
    
        [p1,~,~]              = size(ZPm_plot);
        colorsX               = lines(p1^2);
        color_idx             = 1; 
        legend_entries_real_Z = cell(1, p1);
        legend_entries_imag_Z = cell(1, p1);
        legend_entries_S      = cell(1, p1);
    
        for ind1 = 1:p1
            for ind2 = ind1
                current_color = colorsX(color_idx, :);
                legend_entries_real_Z{color_idx} = ['$',mode,' \mathrm{Re}(Z_{' num2str(ind1) ',' num2str(ind2) '})$'];
                legend_entries_imag_Z{color_idx} = ['$',mode,' \mathrm{Im}(Z_{' num2str(ind1) ',' num2str(ind2) '})$'];
                legend_entries_S{color_idx}      = ['$',mode,' |S_{' num2str(ind1) ',' num2str(ind2) '}|$ (dB)'];
                
                subplot(1,3,1); hold on;
                plot(freq_range/1e+6, real(squeeze(ZPm_plot(ind1,ind2,:))), 'Color', current_color,'LineStyle','-','LineWidth',2);
                subplot(1,3,2); hold on;
                plot(freq_range/1e+6, imag(squeeze(ZPm_plot(ind1,ind2,:))), 'Color', current_color,'LineStyle','-','LineWidth',2);
                subplot(1,3,3); hold on;
                plot(freq_range/1e+6, 20*log10(abs(squeeze(SPm_plot(ind1,ind2,:)))), 'Color', current_color,'LineStyle','-','LineWidth',2);
                
                color_idx = color_idx + 1; 
            end
        end
    
        subplot(1,3,1);
        hx = xline(freq_range(freq_idx)/1e+6); 
        hy = yline(50);
        axis tight;
        ylim([-5 60]);
        xlabel('Frequency (MHz)', 'Interpreter', 'latex', 'FontSize', 14);
        legend(legend_entries_real_Z, 'Interpreter', 'latex', 'FontSize', 14);
        set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
        title('$\mathrm{Re}(Z)$', 'Interpreter', 'latex', 'FontSize', 18);
    
        subplot(1,3,2);
        hx = xline(freq_range(freq_idx)/1e+6); 
        hy = yline(0);
        axis tight;
        ylim([-30 30]);
        xlabel('Frequency (MHz)', 'Interpreter', 'latex', 'FontSize', 14);
        legend(legend_entries_imag_Z, 'Interpreter', 'latex', 'FontSize', 14);
        set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
        title('$\mathrm{Im}(Z)$', 'Interpreter', 'latex', 'FontSize', 18);
    
        subplot(1,3,3);
        hx = xline(freq_range(freq_idx)/1e+6);
        ylim([-80 0]);
        xlabel('Frequency (MHz)', 'Interpreter', 'latex', 'FontSize', 14);
        legend(legend_entries_S, 'Interpreter', 'latex', 'FontSize', 14);
        set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
        title('$\mathrm{dB}(S)$', 'Interpreter', 'latex', 'FontSize', 18);
        
        hold off;

    end

end
function set_black_theme(fig)

    if nargin < 1 || isempty(fig)
        fig = gcf;
    end

    % Set figure background
    set(fig, 'Color', 'k');

    % Loop over all axes in the figure
    ax_all = findall(fig, 'Type', 'axes');
    for ax = ax_all'
        % Axes background and ticks
        set(ax, ...
            'Color', 'k', ...
            'XColor', 'w', ...
            'YColor', 'w', ...
            'ZColor', 'w');  % include Z for 3D plots

        % Labels and title
        set(get(ax, 'XLabel'), 'Color', 'w');
        set(get(ax, 'YLabel'), 'Color', 'w');
        set(get(ax, 'ZLabel'), 'Color', 'w');
        set(get(ax, 'Title'),  'Color', 'w');
    end

    % Legends
    lgd_all = findall(fig, 'Type', 'legend');
    for lgd = lgd_all'
        set(lgd, ...
            'TextColor', 'w', ...
            'Color', 'k', ...
            'EdgeColor', 'w');
    end

    % Colorbars
    cb_all = findall(fig, 'Type', 'colorbar');
    for cb = cb_all'
        set(cb, ...
            'Color', 'w', ...
            'XColor', 'w', ...
            'YColor', 'w');
    end

    % Lines (xline, yline, plot objects, etc.)
    line_objs = findall(fig, 'Type', 'line');
    for ln = line_objs'
        if strcmp(get(ln, 'Color'), [0 0 0])  % change only if black
            set(ln, 'Color', 'w');
        end
    end

    % Text objects (annotations, tick labels, etc.)
    txt_objs = findall(fig, 'Type', 'text');
    for txt = txt_objs'
        set(txt, 'Color', 'w');
    end
end

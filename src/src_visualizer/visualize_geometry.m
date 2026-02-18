function visualize_geometry(MREDM)

    fs = 15;

    dim = MREDM.dimensions;
    WIE = MREDM.WIE;
    SIE = MREDM.SIE;
    
    % ==========================================================
    % Body
    % ==========================================================
    xd = dim.r(:,:,:,1);
    yd = dim.r(:,:,:,2);
    zd = dim.r(:,:,:,3);
    id = dim.idxS;
    P  = [xd(id) yd(id) zd(id)];
    k  = boundary(P,1);        

    % ===============================
    % Laplacian smoothing (geometry)
    % ===============================
    nIter = 2;        % 10–30 usually good
    lambda = 0.5;      % smoothing strength (0–1)
    N = size(P,1);
    A = sparse(N,N);
    for i=1:size(k,1)
        tri = k(i,:);
        A(tri(1),tri(2)) = 1; A(tri(2),tri(1)) = 1;
        A(tri(1),tri(3)) = 1; A(tri(3),tri(1)) = 1;
        A(tri(2),tri(3)) = 1; A(tri(3),tri(2)) = 1;
    end
    deg = sum(A,2);
    % iterative Laplacian smoothing
    Ps = P;
    for it = 1:nIter
        neighbor_avg = (A*Ps)./deg;
        Ps = Ps + lambda*(neighbor_avg - Ps);
    end
    trisurf(k,Ps(:,1),Ps(:,2),Ps(:,3),'LineStyle','none','FaceColor',[254, 227, 212]./255,'FaceAlpha',1); hold on;
    
    % ==========================================================
    % Wire
    % ==========================================================
    if isfield(WIE,'coil')
        if ~isempty(WIE.coil)
            xw = [WIE.coil.F_point(:, 1)'; WIE.coil.S_point(:, 1)'; NaN(1, size(WIE.coil.F_point, 1))];
            yw = [WIE.coil.F_point(:, 2)'; WIE.coil.S_point(:, 2)'; NaN(1, size(WIE.coil.F_point, 1))];
            zw = [WIE.coil.F_point(:, 3)'; WIE.coil.S_point(:, 3)'; NaN(1, size(WIE.coil.F_point, 1))];
            xw = xw(:);
            yw = yw(:);
            zw = zw(:);
            plot3(xw,yw,zw,'Color',[184, 115, 51]/255,'LineStyle','-','LineWidth',2); hold on;

            nPorts = length(WIE.coil.port);
            pos = zeros(nPorts,3);

            for i = 1:nPorts
                pos(i,:) = WIE.coil.F_point(WIE.coil.port(i).t(1),:);
                c = [0 0 0];
                if strcmp(WIE.coil.port(i).type,'port')
                    c = [1 0 0];
                elseif strcmp(WIE.coil.port(i).type,'element')
                    if isfield(WIE.coil.port(i),'load') && ~isempty(WIE.coil.port(i).load)
                        switch WIE.coil.port(i).load
                            case 'capacitor'
                                c = [1 1 0];          % yellow
                            case 'inductor'
                                c = [0 0.447 0.741];  % blue (MATLAB default blue nicer)
                            case 'mutual_inductor'
                                c = [0.494 0.184 0.556]; % purple
                            case 'resistor'
                                c = [0 0 0];          % black
                            otherwise
                                c = [0.5 0.5 0.5];    % unknown → gray
                        end
                    end
                end
                scatter3(pos(i,1),pos(i,2),pos(i,3), 120,'o', 'MarkerFaceColor',c,'MarkerEdgeColor',c,'LineWidth',1); hold on;
            end
        end
    end
    
    % ==========================================================
    % Surface Coil
    % ==========================================================
    if isfield(SIE,'coil')
        if ~isempty(SIE.coil)
            TR1 = triangulation(SIE.coil.elem(1:3,:)',SIE.coil.node');
            trisurf(TR1,'FaceColor',[184, 115, 51]/255,'FaceAlpha',1,'LineStyle','-');  hold on;    

            nPorts = length(SIE.coil.port);
            pos = zeros(nPorts,3);

            for i = 1:nPorts
                pos(i,:) = SIE.coil.node(:, SIE.coil.edge(1, SIE.coil.index == SIE.coil.port(i).t(1)));
                c = [0 0 0];
                if strcmp(SIE.coil.port(i).type,'port')
                    c = [1 0 0];
                elseif strcmp(SIE.coil.port(i).type,'element')
                    if isfield(SIE.coil.port(i),'load') && ~isempty(SIE.coil.port(i).load)
                        switch SIE.coil.port(i).load
                            case 'capacitor'
                                c = [1 1 0];          % yellow
                            case 'inductor'
                                c = [0 0.447 0.741];  % blue (MATLAB default blue nicer)
                            case 'mutual_inductor'
                                c = [0.494 0.184 0.556]; % purple
                            case 'resistor'
                                c = [0 0 0];          % black
                            otherwise
                                c = [0.5 0.5 0.5];    % unknown → gray
                        end
                    end
                end
                scatter3(pos(i,1),pos(i,2),pos(i,3), 120,'o', 'MarkerFaceColor',c,'MarkerEdgeColor',c,'LineWidth',1); hold on;
            end

        end
    end

    % ==========================================================
    % Surface Shield
    % ==========================================================
    if isfield(SIE,'shield')
        if ~isempty(SIE.shield)
            TR2 = triangulation(SIE.shield.elem(1:3,:)',SIE.shield.node');
            trisurf(TR2,'FaceColor',[184, 115, 51]/255,'FaceAlpha',0.3,'LineStyle','none');  hold on;

            nPorts = length(SIE.shield.port);
            pos = zeros(nPorts,3);

            for i = 1:nPorts
                pos(i,:) = SIE.shield.node(:, SIE.shield.edge(1, SIE.shield.index == SIE.shield.port(i).t(1)));
                c = [0 0 0];
                if strcmp(SIE.shield.port(i).type,'port')
                    c = [1 0 0];
                elseif strcmp(SIE.shield.port(i).type,'element')
                    if isfield(SIE.shield.port(i),'load') && ~isempty(SIE.shield.port(i).load)
                        switch SIE.shield.port(i).load
                            case 'capacitor'
                                c = [1 1 0];          % yellow
                            case 'inductor'
                                c = [0 0.447 0.741];  % blue (MATLAB default blue nicer)
                            case 'mutual_inductor'
                                c = [0.494 0.184 0.556]; % purple
                            case 'resistor'
                                c = [0 0 0];          % black
                            otherwise
                                c = [0.5 0.5 0.5];    % unknown → gray
                        end
                    end
                end
                scatter3(pos(i,1),pos(i,2),pos(i,3), 120,'o', 'MarkerFaceColor',c,'MarkerEdgeColor',c,'LineWidth',1); hold on;
            end
        end
    end

    % ==========================================================
    % Dipole Basis Support
    % ==========================================================
    if isfield(dim,'r_basis')
        xb = dim.r_basis(:,:,:,1);
        yb = dim.r_basis(:,:,:,2);
        zb = dim.r_basis(:,:,:,3);
        ib = dim.basis_idxS;
        scatter3(xb(ib),yb(ib),zb(ib),20,'sq','MarkerEdgeColor',[255, 165, 0]/255,'MarkerFaceColor',[255, 165, 0]/255,'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);  hold on;
    end
    
    xlabel('x (m)','interpreter','latex','Fontsize',fs);
    ylabel('y (m)','interpreter','latex','Fontsize',fs);
    zlabel('z (m)','interpreter','latex','Fontsize',fs);
    
    ax = gca;
    ax.FontSize = fs;
    axis equal;
    axis on;
    grid on;
    lighting gouraud    
    camlight headlight   
    material dull        

    hold off;
    
end
    
% pos = zeros(length(COIL.port),3);
% s = zeros(length(COIL.port),1);
% for i = 1:length(COIL.port)
%     pos(i,:) = COIL.node(:,COIL.edge(1,COIL.index == COIL.port(i).t(1)));
%     s(i) = i;
% end
% s = string(s);
% text(pos(:,1),pos(:,2),pos(:,3),s,'Interpreter','none','HorizontalAlignment','center','Color','k','FontSize',fs+5,'FontWeight','bold'); hold on;

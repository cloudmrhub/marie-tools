function [Z,Z_loss] = assembly_le(Z,Z_loss,coil,port,emc)
    
    index = coil.index;
    node = coil.node;
    edge = coil.edge;
    
    Nterms = length(port.t);

    if strcmp(port.load,'mutual_inductor')
        N2terms = length(port.t_mutual);
    end
    
    %% select load type
    if strcmp(port.load,'resistor')
        Zx_val = port.value;
        R_loss = port.value;  
    elseif strcmp(port.load,'inductor')
        Zx_val = 1i * emc.omega * port.value;
        R_loss = (emc.omega * port.value) / port.Q;
        Zx_val = Zx_val+R_loss;
    elseif strcmp(port.load,'capacitor')
        Zx_val = 1 / (1i * emc.omega * port.value);
        R_loss = 1 / (emc.omega * port.value * port.Q);
        Zx_val = Zx_val+R_loss;
    elseif strcmp(port.load,'mutual_inductor')
        Zx_val = 1i * emc.omega * port.value;
        Zx_val_M = 1i * emc.omega * port.value_mutual;
        R_loss = (emc.omega * port.value) / port.Q;
        Zx_val = Zx_val+R_loss;
    end
	
    %% update Method of Moments matrix
    L_edges = zeros(Nterms,1);
    for i = 1:Nterms
        index_edge = port.t(i);
        phys_edge = find(index == index_edge);
        point1 = edge(1,phys_edge);
        point2 = edge(2,phys_edge);
        rpoint1 = node(:, point1);
        rpoint2 = node(:, point2);
        L_edges(i) = norm(rpoint1-rpoint2);
    end

    if strcmp(port.load,'mutual_inductor')
        L2_edges = zeros(N2terms,1);
        for i = 1:N2terms
            index_edge = port.t_mutual(i);
            phys_edge = find(index == index_edge);
    	    point1 = edge(1,phys_edge);
            point2 = edge(2,phys_edge);
            rpoint1 = node(:, point1);
            rpoint2 = node(:, point2);
            L2_edges(i) = norm(rpoint1-rpoint2);     
        end
    end

    for i = 1:Nterms
        index_i = port.t(i);
        for j = 1:Nterms
            index_j = port.t(j);
            Z_update = Zx_val * L_edges(i) * L_edges(j);
            Z(index_i, index_j) = Z(index_i, index_j) + Z_update;
            Z_loss(index_i, index_j) = Z_loss(index_i, index_j) + R_loss *  L_edges(i) * L_edges(j);
        end
    end

    if strcmp(port.load,'mutual_inductor')
        for i = 1:Nterms
            index_i = port.t(i);
            for j = 1:N2terms
                index_j = port.t_mutual(j);
                Z_update = Zx_val_M * L_edges(i) * L2_edges(j);
                Z(index_i, index_j) = Z(index_i, index_j) + Z_update;
            end
        end
    end    
end
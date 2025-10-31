function [V] = excitation_coil(coil,port)

    Vmult   = port.voltage;

    index = coil.index;
    node  = coil.node;
    edge  = coil.edge;

	Nsie = max(index);
	V  = zeros(Nsie,1);
	Nfeed_terms = length(port.t);
	
    for i = 1:Nfeed_terms    
        index_edge        = port.t(i);
        real_edge         = find(index == index_edge);  
        point1            = edge(1,real_edge);
        point2            = edge(2,real_edge);   
        rpoint1           = node(:,point1);
        rpoint2           = node(:,point2);     
        L_edge            = norm(rpoint2-rpoint1);
        V(index_edge,1)   = -L_edge*Vmult; 
    end
	
end
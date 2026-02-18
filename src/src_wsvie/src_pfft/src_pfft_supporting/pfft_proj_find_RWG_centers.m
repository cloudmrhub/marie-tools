function [RWG_cntr] = pfft_proj_find_RWG_centers(coil)

    index = coil.index;
    edge  = coil.edge;
    node  = coil.node;

    N_sie = max(index);

    RWG_cntr = zeros(N_sie, 3); 

    for i_sie = 1:N_sie
        
        phys_idx = find(index == i_sie);
        node1 = edge(1, phys_idx);
        node2 = edge(2, phys_idx);
        r1   = node(:,node1);
        r2   = node(:,node2);
        r_c  = (r1(:) + r2(:)) / 2;
        RWG_cntr(i_sie, :) = r_c.';
    
    end
    
end

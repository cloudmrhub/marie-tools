function r_rwg = get_rwg_vertices (coil, i_sie)
    % function extracts vertices of the RWG basis function, that corresponds to
    % a given i_sie

    % find physical edge for given dof(i_sie)
    phys_edge = find(coil.index == i_sie);

    % find tirangles that share this physical edge
    [idx_tr1,tr_1] = find(coil.etod == phys_edge);
    [idx_tr2,tr_2] = find(coil.etod == -phys_edge);

    % get coordinates of oposit vertices
    rv_p = coil.node(:,coil.elem(idx_tr1,tr_1)); 
    rv_n = coil.node(:,coil.elem(idx_tr2,tr_2));

    % get local index of vertices, that form common edge
    edge_node = setdiff([1,2,3], idx_tr1);

    % get coordinates of vertices, that form common edge
    r23 = coil.node(:,coil.elem(edge_node,tr_1));

    % stack coordinates in a single vector
    r_rwg = [rv_p; rv_n; r23(:)];

end
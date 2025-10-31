function [SHIELD] = geo_shield(filename) 

    port = repmat(struct('t', [], 'type',[], 'load', [], 'value', [], 't_mutual', [], 'value_mutual', [], 'voltage', []), 0, 1);
    
    [node,~,e,elem,~] = Mesh_Parse(filename);
    elem = Mesh_Permute(elem,e);

    [edge,etod,index,~,index_elem] = Mesh_PreProc(e,elem,port);
    [Ct,Ln,Pn] = Mesh_CLP(node,elem);

    SHIELD = struct('index',index,'etod',etod,'node',node,'edge',edge,'elem',elem,'index_elem',index_elem,'Ct',Ct,'Ln',Ln,'Pn',Pn);
      
    N_sie = max(index);
    r_v = zeros(12, N_sie);
    for i = 1:N_sie
        r_v(:,i) = get_rwg_vertices(SHIELD, i);
    end         
    r_v = r_v.';
    SHIELD.rp = r_v(:,1:3).';
    SHIELD.rn = r_v(:,4:6).';
    SHIELD.r2 = r_v(:,7:9).';
    SHIELD.r3 = r_v(:,10:12).';
end
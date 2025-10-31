function [SCOIL] = geo_scoil(filename,LC) 

    Nports = length(LC);
    portstruct = struct('t', [], 'type',[], 'load', [], 'value', [], 't_mutual', [], 'value_mutual', [], 'voltage', []);
    port_F = repmat(portstruct, Nports, 1);
    
    for port_num = 1:Nports
    
        port_F(port_num).type = LC(port_num).type;
        port_F(port_num).load = LC(port_num).load;
        port_F(port_num).value = LC(port_num).value;
        port_F(port_num).voltage = LC(port_num).excitation.voltage;
        
        if strcmp(LC(port_num).load,'mutual_inductor')
            port_F(port_num).t_mutual = LC(port_num).cross_talk(1);
			port_F(port_num).value_mutual = LC(port_num).cross_talk(2);
        end
        
    end
    
    [node,~,e,elem,~] = Mesh_Parse(filename);
    elem = Mesh_Permute(elem,e);

    [edge,etod,index,port,index_elem] = Mesh_PreProc(e,elem,port_F);

    for port_num = 1:Nports
        if strcmp(port(port_num).load,'mutual_inductor')
			port(port_num).t_mutual = port(port(port_num).t_mutual).t;
        end 
    end

    [Ct,Ln,Pn] = Mesh_CLP(node,elem);

    SCOIL = struct('index',index,'etod',etod,'node',node,'edge',edge,'elem',elem,'index_elem',index_elem,'Ct',Ct,'Ln',Ln,'Pn',Pn,'port',port);
        
    N_sie = max(index);
    r_v = zeros(12, N_sie);
    for i = 1:N_sie
        r_v(:,i) = get_rwg_vertices(SCOIL, i);
    end         
    r_v = r_v.';
    SCOIL.rp = r_v(:,1:3).';
    SCOIL.rn = r_v(:,4:6).';
    SCOIL.r2 = r_v(:,7:9).';
    SCOIL.r3 = r_v(:,10:12).';
end
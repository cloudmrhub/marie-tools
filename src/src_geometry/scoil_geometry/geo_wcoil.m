function [WCOIL] = geo_wcoil(filename,LC) 

    if strcmp(filename(end-2:end),'wsg')
        S = load(filename,'-ascii');
        portentry        = find(S(:,7)>0);
        loop_start       = find(squeeze(S(:,8))==1);
        loop_end         = find(squeeze(S(:,8))==2);
        F_point          = S(:,1:3);
        S_point          = S(:,4:6);
        T_point          = S(2:end,4:6);
    elseif strcmp(filename(end-2:end),'msh')
       [portentry,...
        loop_start,...
        loop_end,...
        F_point,...
        S_point,...
        T_point] = Mesh_Wire(filename);
    end

    Nports = length(LC);
    portstruct = struct('t', [], 'type',[], 'load', [], 'value', [], 't_mutual', [], 'value_mutual', [], 'voltage', []);
    port = repmat(portstruct, Nports, 1);
    for port_num = 1:Nports
        port(port_num).t    = portentry(port_num);
        port(port_num).type = LC(port_num).type;
        port(port_num).load = LC(port_num).load;
        port(port_num).value = LC(port_num).value;
        port(port_num).Q = LC(port_num).Q;
        port(port_num).voltage = LC(port_num).excitation.voltage;
        if strcmp(LC(port_num).load,'mutual_inductor')
            port(port_num).t_mutual = LC(port_num).cross_talk(1);
			port(port_num).value_mutual = LC(port_num).cross_talk(2);
        end        
    end

    for port_num = 1:Nports
        if strcmp(port(port_num).load,'mutual_inductor')
			port(port_num).t_mutual = port(port(port_num).t_mutual).t;
        end 
    end

    WCOIL = struct('F_point',F_point,'S_point',S_point,'T_point',T_point,'loop_start',loop_start,'loop_end',loop_end,'port',port);
    WCOIL = ProcessLoops(WCOIL);

end

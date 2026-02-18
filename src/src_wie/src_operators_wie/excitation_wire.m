function [F] = excitation_wire(Wcoil,N_wie,N_ports,N_loads)

	F   = zeros(N_wie,N_ports);
    port_counter = 0;
    for i = 1:N_loads
        if strcmp(Wcoil.port(i).type,'port')
            V   = Wcoil.port(i).voltage;
            port_counter = port_counter + 1;
            s1 = Wcoil.port(i).t;
            s2 = Wcoil.port(i).t+1;
            F(s1,port_counter)   = -0.5*V;
            F(s2,port_counter)   = -0.5*V;
        end
    end
   
    if isempty(Wcoil.loop_end)
        F(1:end-1,:) = F(2:end,:);
        F = F(1:end-1,:);
    end
	
end

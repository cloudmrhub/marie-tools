function [Z,Z_loss] = assembly_wle(Z,Z_loss,port,emc)
    
    %% Select load type
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
	
    %% Update Method of Moments matrix
    s1 = port.t;
    s2 = port.t+1;
    Z(s1,s1) = Z(s1,s1) + 0.5 * Zx_val;
    Z(s2,s2) = Z(s2,s2) + 0.5 * Zx_val;

    %% Update resistance-only matrix
    Z_loss(s1, s1) = Z_loss(s1, s1) + 0.5 * R_loss;
    Z_loss(s2, s2) = Z_loss(s2, s2) + 0.5 * R_loss;
    
    if strcmp(port.load,'mutual_inductor')
        s11 = port.t;
        s21 = port.t_mutual;
        s12 = port.t+1;
        s22 = port.t_mutual+1;
        
        Z(s11,s21) = Z(s11,s21) + 0.5*Zx_val_M;
        Z(s12,s22) = Z(s12,s22) + 0.5*Zx_val_M;
    end    

end
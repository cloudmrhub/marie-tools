function[MREDM] = co_simulation_store_elements(MREDM,separator_counter,shield_file,wire_file,coil_file,RLC_org_h,RLC_tmd_h,RLC_org_w,RLC_tmd_w,RLC_org_s,RLC_tmd_s,final_elements) 

    if ~isempty(MREDM.WIE.RLC) && isempty(MREDM.SIE.RLC) && isempty(MREDM.SIE.shield_RLC)
        MREDM.dimensions.N_ports = geo_scoil_print_lumped_elements(wire_file,RLC_org_w,RLC_tmd_w,final_elements);    

    elseif isempty(MREDM.WIE.RLC) && ~isempty(MREDM.SIE.RLC) && isempty(MREDM.SIE.shield_RLC)
        MREDM.dimensions.M_ports = geo_scoil_print_lumped_elements(coil_file,RLC_org_s,RLC_tmd_s,final_elements);
    
    elseif ~isempty(MREDM.WIE.RLC) && ~isempty(MREDM.SIE.RLC) && isempty(MREDM.SIE.shield_RLC)
        wire_elements = final_elements(1:separator_counter);
        coil_elements = final_elements(separator_counter+1:end);
        MREDM.dimensions.M_ports = geo_scoil_print_lumped_elements(wire_file,RLC_org_w,RLC_tmd_w,wire_elements);
        MREDM.dimensions.N_ports = geo_scoil_print_lumped_elements(coil_file,RLC_org_s,RLC_tmd_s,coil_elements);

    elseif ~isempty(MREDM.WIE.RLC) && isempty(MREDM.SIE.RLC) && ~isempty(MREDM.SIE.shield_RLC)
        shield_elements = final_elements(1:separator_counter);
        wire_elements   = final_elements(separator_counter+1:end);
        MREDM.dimensions.M_ports        = geo_scoil_print_lumped_elements(wire_file,RLC_org_w,RLC_tmd_w,wire_elements);
        MREDM.dimensions.N_ports_shield = geo_scoil_print_lumped_elements(shield_file,RLC_org_h,RLC_tmd_h,shield_elements);

    elseif isempty(MREDM.WIE.RLC) && ~isempty(MREDM.SIE.RLC) && ~isempty(MREDM.SIE.shield_RLC)
        shield_elements = final_elements(1:separator_counter);
        coil_elements   = final_elements(separator_counter+1:end);
        MREDM.dimensions.N_ports        = geo_scoil_print_lumped_elements(coil_file,RLC_org_s,RLC_tmd_s,coil_elements);
        MREDM.dimensions.N_ports_shield = geo_scoil_print_lumped_elements(shield_file,RLC_org_h,RLC_tmd_h,shield_elements);
    
    elseif ~isempty(MREDM.WIE.RLC) && ~isempty(MREDM.SIE.RLC) && ~isempty(MREDM.SIE.shield_RLC)
        shield_elements = final_elements(1:separator_counter.first);
        wire_elements   = final_elements(separator_counter.first+1:separator_counter.first+separator_counter.second);
        coil_elements   = final_elements(separator_counter.first+separator_counter.second+1:end);
        MREDM.dimensions.M_ports        = geo_scoil_print_lumped_elements(wire_file,RLC_org_w,RLC_tmd_w,wire_elements);
        MREDM.dimensions.N_ports        = geo_scoil_print_lumped_elements(coil_file,RLC_org_s,RLC_tmd_s,coil_elements);
        MREDM.dimensions.N_ports_shield = geo_scoil_print_lumped_elements(shield_file,RLC_org_h,RLC_tmd_h,shield_elements);
        
    end
end
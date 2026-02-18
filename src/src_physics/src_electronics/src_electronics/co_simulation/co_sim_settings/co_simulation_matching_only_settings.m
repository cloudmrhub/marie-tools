function[matching_lumped_elements, port_order, loads_matching] = co_simulation_matching_only_settings(RLC_org)

    % Set up optimization arrays
    matching_element_counter = 0;
    port_counter = 0;
    for i = 1:length(RLC_org)
        if strcmp(RLC_org(i).type,'port')
            port_counter = port_counter+1;
            matching_element_counter = matching_element_counter + length(convertCharsToStrings(RLC_org(i).load));
        end
    end

    loads_matching = strings(matching_element_counter,1);
    port_order     = zeros(matching_element_counter,1);

    c1 = 0;
    matching_lumped_elements = [];
    for port_num = 1:length(RLC_org)
        port_order_counter = 1;
        if strcmp(RLC_org(port_num).type,'port')
            for i = 1:length(convertCharsToStrings(RLC_org(port_num).load))
                matching_lumped_elements(c1+1)             = RLC_org(port_num).value(i);
                load_tmp                                   = convertCharsToStrings(RLC_org(port_num).load);
                loads_matching(c1+1)                       = load_tmp(i);
                port_order(c1+1)                           = port_order_counter;
                c1 = c1+1;
                port_order_counter = port_order_counter+1;
            end
        end
    end

end
function[match_split,...
         mask_detailed_split] = co_simulation_split_match_settings(RLC_org)
    
    c1 = 1;
    for port_num = 1:length(RLC_org)
        if strcmp(RLC_org(port_num).type,'port')
            if strcmp(RLC_org(port_num).excitation.TxRx,'Rx')
                match_split(c1) = 1;
            elseif strcmp(RLC_org(port_num).excitation.TxRx,'Tx')
                match_split(c1) = 2;
            elseif strcmp(RLC_org(port_num).excitation.TxRx,'TxRx')
                match_split(c1) = 3;
            end
            c1 = c1 + 1;
        end
    end

    c3 = 1;
    for port_num = 1:length(RLC_org)
        if strcmp(RLC_org(port_num).type,'port')
            for i = 1:numel(RLC_org(1).load)
                if strcmp(RLC_org(port_num).excitation.TxRx,'Rx')
                    mask_detailed_split(c3) = 1;
                elseif strcmp(RLC_org(port_num).excitation.TxRx,'Tx')
                    mask_detailed_split(c3) = 2;
                elseif strcmp(RLC_org(port_num).excitation.TxRx,'TxRx')
                    mask_detailed_split(c3) = 3;
                end
                c3 = c3 + 1;
            end
        end
    end

    match_split         = match_split';
    mask_detailed_split = mask_detailed_split';

end
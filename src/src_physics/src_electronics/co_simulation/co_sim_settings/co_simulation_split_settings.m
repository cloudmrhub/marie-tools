function[match_split,...
         tune_split,...
         mask_split,...
         matching_mask_split,...
         fixed_matching_mask_split,...
         tuning_mask_split,...
         mutual_mask_split,...
         mask_detailed_split,...
         fixed_mask_detailed_split] = co_simulation_split_settings(RLC_tmd,...
                                                                   RLC_org)

    % Get variable elements
    tuning_elements   = 0;
    matching_elements = 0;
    mutual_elements   = 0;
    for port_num = 1:length(RLC_org)
        if (RLC_tmd(port_num).optim.boolean == 1 && strcmp(RLC_org(port_num).type,'port'))  || (RLC_tmd(port_num).optim.boolean == 0 && strcmp(RLC_org(port_num).type,'port'))
            matching_elements = matching_elements+1;
        elseif RLC_tmd(port_num).optim.boolean == 1 && strcmp(RLC_org(port_num).type,'element')
            tuning_elements = tuning_elements+1;
            if strcmp(RLC_org(port_num).load,'mutual_inductor')
                mutual_elements = mutual_elements+1/2;
            end
        end
    end
    mutual_elements           = round(mutual_elements);
    all_elements              = matching_elements+tuning_elements+mutual_elements;
    matching_mask_split       = zeros(all_elements,1);
    fixed_matching_mask_split = zeros(all_elements,1);
    tuning_mask_split         = zeros(all_elements,1);
    mutual_mask_split         = zeros(all_elements,1);
    c1 = 0;
    c2 = 0;
    c3 = 0;
    c4 = 0;
    mutual_mod = 1;
    for port_num = 1:length(RLC_org)
        if RLC_tmd(port_num).optim.boolean == 1  || (RLC_tmd(port_num).optim.boolean == 0 && strcmp(RLC_org(port_num).type,'port'))
            if strcmp(RLC_org(port_num).type,'port') && RLC_tmd(port_num).optim.boolean == 1
                c1 = c1+1;
                if strcmp(RLC_tmd(port_num).excitation.TxRx,'Rx')
                    matching_mask_split(c1) = 1;
                elseif strcmp(RLC_tmd(port_num).excitation.TxRx,'Tx')
                    matching_mask_split(c1) = 2;
                elseif strcmp(RLC_tmd(port_num).excitation.TxRx,'TxRx')
                    matching_mask_split(c1) = 3;
                end
            elseif strcmp(RLC_org(port_num).type,'port') && RLC_tmd(port_num).optim.boolean == 0
                c4 = c4+1;
                if strcmp(RLC_tmd(port_num).excitation.TxRx,'Rx')
                    fixed_matching_mask_split(c4) = 1;
                elseif strcmp(RLC_tmd(port_num).excitation.TxRx,'Tx')
                    fixed_matching_mask_split(c4) = 2;
                elseif strcmp(RLC_tmd(port_num).excitation.TxRx,'TxRx')
                    fixed_matching_mask_split(c4) = 3;
                end
            elseif strcmp(RLC_org(port_num).type,'element')
                c2 = c2+1;
                if strcmp(RLC_tmd(port_num).excitation.TxRx,'Rx')
                    tuning_mask_split(c2) = 1;
                elseif strcmp(RLC_tmd(port_num).excitation.TxRx,'Tx')
                    tuning_mask_split(c2) = 2;
                elseif strcmp(RLC_tmd(port_num).excitation.TxRx,'TxRx')
                    tuning_mask_split(c2) = 3;
                end
                if strcmp(RLC_org(port_num).load,'mutual_inductor')
                    if mod(mutual_mod,2)
                        c3 = c3+1;
                        if strcmp(RLC_tmd(port_num).excitation.TxRx,'Rx')
                            mutual_mask_split(c3) = 1;
                        elseif strcmp(RLC_tmd(port_num).excitation.TxRx,'Tx')
                            mutual_mask_split(c3) = 2;
                        elseif strcmp(RLC_tmd(port_num).excitation.TxRx,'TxRx')
                            mutual_mask_split(c3) = 3;
                        end
                    end
                    mutual_mod = mutual_mod + 1;
                end
            end
        end
    end 
    
    c1 = 1;
    for port_num = 1:length(RLC_org)
        if (RLC_tmd(port_num).optim.boolean == 1 && strcmp(RLC_org(port_num).type,'port'))  || (RLC_tmd(port_num).optim.boolean == 0 && strcmp(RLC_org(port_num).type,'port'))
            if strcmp(RLC_tmd(port_num).excitation.TxRx,'Rx')
                match_split(c1) = 1;
            elseif strcmp(RLC_tmd(port_num).excitation.TxRx,'Tx')
                match_split(c1) = 2;
            elseif strcmp(RLC_tmd(port_num).excitation.TxRx,'TxRx')
                match_split(c1) = 3;
            end
            c1 = c1 + 1;
        end
    end
    
    c2 = 1;
    for port_num = 1:length(RLC_org)
        if RLC_tmd(port_num).optim.boolean == 1 && strcmp(RLC_org(port_num).type,'element')
            if strcmp(RLC_tmd(port_num).excitation.TxRx,'Rx')
                tune_split(c2) = 1;
            elseif strcmp(RLC_tmd(port_num).excitation.TxRx,'Tx')
                tune_split(c2) = 2;
            elseif strcmp(RLC_tmd(port_num).excitation.TxRx,'TxRx')
                tune_split(c2) = 3;
            end
            c2 = c2 + 1;
        end
    end

    c3 = 1;
    for port_num = 1:length(RLC_org)
        if (RLC_tmd(port_num).optim.boolean == 1) || (RLC_tmd(port_num).optim.boolean == 0 && strcmp(RLC_org(port_num).type,'port'))
            if strcmp(RLC_org(port_num).type,'port')
                for i = 1:numel(RLC_tmd(port_num).load)
                    if strcmp(RLC_tmd(port_num).excitation.TxRx,'Rx')
                        mask_detailed_split(c3) = 1;
                        if RLC_tmd(port_num).optim.boolean == 0
                            fixed_mask_detailed_split(c3) = false;
                        elseif RLC_tmd(port_num).optim.boolean == 1
                            fixed_mask_detailed_split(c3) = true;
                        end
                    elseif strcmp(RLC_tmd(port_num).excitation.TxRx,'Tx')
                        mask_detailed_split(c3) = 2;
                        if RLC_tmd(port_num).optim.boolean == 0
                            fixed_mask_detailed_split(c3) = false;
                        elseif RLC_tmd(port_num).optim.boolean == 1
                            fixed_mask_detailed_split(c3) = true;
                        end
                    elseif strcmp(RLC_tmd(port_num).excitation.TxRx,'TxRx')
                        mask_detailed_split(c3) = 3;
                        if RLC_tmd(port_num).optim.boolean == 0
                            fixed_mask_detailed_split(c3) = false;
                        elseif RLC_tmd(port_num).optim.boolean == 1
                            fixed_mask_detailed_split(c3) = true;
                        end
                    end
                    c3 = c3 + 1;
                end
            else
                if strcmp(RLC_tmd(port_num).excitation.TxRx,'Rx')
                    mask_detailed_split(c3) = 1;
                    fixed_mask_detailed_split(c3) = true;
                elseif strcmp(RLC_tmd(port_num).excitation.TxRx,'Tx')
                    mask_detailed_split(c3) = 2;
                    fixed_mask_detailed_split(c3) = true;
                elseif strcmp(RLC_tmd(port_num).excitation.TxRx,'TxRx')
                    mask_detailed_split(c3) = 3;
                    fixed_mask_detailed_split(c3) = true;
                end
                c3 = c3 + 1;
            end
        end
    end

    c4 = 1;
    for port_num = 1:length(RLC_org)
        if RLC_tmd(port_num).optim.boolean == 1 || (RLC_tmd(port_num).optim.boolean == 0 && strcmp(RLC_org(port_num).type,'port'))
            if strcmp(RLC_tmd(port_num).excitation.TxRx,'Rx')
                mask_split(c4) = 1;
            elseif strcmp(RLC_tmd(port_num).excitation.TxRx,'Tx')
                mask_split(c4) = 2;
            elseif strcmp(RLC_tmd(port_num).excitation.TxRx,'TxRx')
                mask_split(c4) = 3;
            end
            c4 = c4 + 1;
        end
    end

    match_split         = match_split';
    mask_split          = mask_split';
    mask_detailed_split = mask_detailed_split';
    fixed_mask_detailed_split = fixed_mask_detailed_split.';

end
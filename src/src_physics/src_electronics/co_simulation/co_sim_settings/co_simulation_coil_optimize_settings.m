function[RLC_org,RLC_tmd,RLC_org_h,RLC_tmd_h,RLC_org_w,RLC_tmd_w,RLC_org_s,RLC_tmd_s,separator_counter,shield_file,wire_file,coil_file] = co_simulation_coil_optimize_settings(MREDM)  
    
    separator_counter = 0;
    wire_file = [];
    coil_file = [];
    shield_file = [];

    % Coil
    if isempty(MREDM.WIE.RLC) && ~isempty(MREDM.SIE.RLC) && isempty(MREDM.SIE.shield_RLC)
        RLC_tmd   = MREDM.SIE.RLC;
        coil_file = fullfile('./data/coils/coil_files',MREDM.inputs.coil);
        RLC_org   = geo_scoil_lumped_elements(coil_file,0);
        RLC_tmd_s = RLC_tmd;
        RLC_org_s = RLC_org;
        RLC_tmd_w = [];
        RLC_org_w = [];
        RLC_tmd_h = [];
        RLC_org_h = [];
    % Wire
    elseif ~isempty(MREDM.WIE.RLC) && isempty(MREDM.SIE.RLC) && isempty(MREDM.SIE.shield_RLC)
        RLC_tmd   = MREDM.WIE.RLC;
        wire_file = fullfile('./data/coils/wire_files',MREDM.inputs.wire);
        RLC_org   = geo_scoil_lumped_elements(wire_file,0);
        RLC_tmd_w = RLC_tmd;
        RLC_org_w = RLC_org;
        RLC_tmd_s = [];
        RLC_org_s = [];
        RLC_tmd_h = [];
        RLC_org_h = [];
    % Shield
    elseif isempty(MREDM.WIE.RLC) && isempty(MREDM.SIE.RLC) && ~isempty(MREDM.SIE.shield_RLC)
        RLC_tmd     = MREDM.SIE.shield_RLC;
        shield_file = fullfile('./data/coils/shield_files',MREDM.inputs.shield);
        RLC_org     = geo_scoil_lumped_elements(shield_file,0);
        RLC_tmd_h   = RLC_tmd;
        RLC_org_h   = RLC_org;
        RLC_tmd_s   = [];
        RLC_org_s   = [];
        RLC_tmd_w   = [];
        RLC_org_w   = [];
    % Wire - Coil
    elseif ~isempty(MREDM.WIE.RLC) && ~isempty(MREDM.SIE.RLC) && isempty(MREDM.SIE.shield_RLC)
        RLC_tmd_s = MREDM.SIE.RLC;
        RLC_tmd_w = MREDM.WIE.RLC;
        coil_file = fullfile('./data/coils/coil_files',MREDM.inputs.coil);
        RLC_org_s = geo_scoil_lumped_elements(coil_file,0);
        wire_file = fullfile('./data/coils/wire_files',MREDM.inputs.wire);
        RLC_org_w = geo_scoil_lumped_elements(wire_file,0);
        nW        = RLC_org_w(end).number;
        nS        = RLC_org_s(end).number;
        RLC_org   = repmat(struct('number',[], 'type',[], 'load',[], 'value',[], 'optim',[], 'cross_talk', [], 'excitation', []), nW+nS,1);
        RLC_tmd   = repmat(struct('number',[], 'type',[], 'load',[], 'value',[], 'optim',[], 'cross_talk', [], 'excitation', []), nW+nS,1);
        % Concatenate
        RLC_org(1:nW) = RLC_org_w(:);
        RLC_tmd(1:nW) = RLC_tmd_w(:);
        RLC_org(nW+1:end) = RLC_org_s(:);
        RLC_tmd(nW+1:end) = RLC_tmd_s(:);
        separator_counter = 0;
        for i = 1:length(RLC_tmd_w)
            if RLC_tmd_w(i).optim.boolean == 1
                separator_counter = separator_counter + length(convertCharsToStrings(RLC_tmd_w(i).load));
            end
        end
        max_entity_org = max(arrayfun(@(x) x.excitation.entity, RLC_org_w));
        max_entity_tmd = max(arrayfun(@(x) x.excitation.entity, RLC_tmd_w));
        max_symmetry = max(arrayfun(@(x) x.optim.symmetry, RLC_tmd_w));
        % Fix ordering
        for i = 1:nS
            if ~isempty(RLC_org_s(i).cross_talk)
                RLC_org(nW+i).cross_talk(2) = RLC_org_s(i).cross_talk(2)+nW;
            end
            if ~isempty(RLC_tmd_s(i).cross_talk)
                RLC_tmd(nW+i).cross_talk(2) = RLC_tmd_s(i).cross_talk(2)+nW;
            end
            RLC_org(nW+i).number            = RLC_org_s(i).number+nW;
            RLC_tmd(nW+i).number            = RLC_tmd_s(i).number+nW;
            RLC_org(nW+i).excitation.entity = RLC_org_s(i).excitation.entity + max_entity_org;
            RLC_tmd(nW+i).excitation.entity = RLC_tmd_s(i).excitation.entity + max_entity_tmd;
            if RLC_tmd(nW+i).optim.boolean == 1
                RLC_tmd(nW+i).optim.symmetry = RLC_tmd_s(i).optim.symmetry + max_symmetry;
            end
        end
        RLC_tmd_h = [];
        RLC_org_h = [];
    % Coil - Shield
    elseif isempty(MREDM.WIE.RLC) && ~isempty(MREDM.SIE.RLC) && ~isempty(MREDM.SIE.shield_RLC)
        RLC_tmd_s = MREDM.SIE.RLC;
        RLC_tmd_h = MREDM.SIE.shield_RLC;
        coil_file = fullfile('./data/coils/coil_files',MREDM.inputs.coil);
        RLC_org_s = geo_scoil_lumped_elements(coil_file,0);
        shield_file = fullfile('./data/coils/shield_files',MREDM.inputs.shield);
        RLC_org_h = geo_scoil_lumped_elements(shield_file,0);
        nH        = RLC_org_h(end).number;
        nS        = RLC_org_s(end).number;
        RLC_org   = repmat(struct('number',[], 'type',[], 'load',[], 'value',[], 'optim',[], 'cross_talk', [], 'excitation', []), nH+nS,1);
        RLC_tmd   = repmat(struct('number',[], 'type',[], 'load',[], 'value',[], 'optim',[], 'cross_talk', [], 'excitation', []), nH+nS,1);
        % Concatenate
        RLC_org(1:nH) = RLC_org_h(:);
        RLC_tmd(1:nH) = RLC_tmd_h(:);
        RLC_org(nH+1:end) = RLC_org_s(:);
        RLC_tmd(nH+1:end) = RLC_tmd_s(:);
        separator_counter = 0;
        for i = 1:length(RLC_tmd_h)
            if RLC_tmd_h(i).optim.boolean == 1
                separator_counter = separator_counter + length(convertCharsToStrings(RLC_tmd_h(i).load));
            end
        end
        max_entity_org = max(arrayfun(@(x) x.excitation.entity, RLC_org_h));
        max_entity_tmd = max(arrayfun(@(x) x.excitation.entity, RLC_tmd_h));
        max_symmetry = max(arrayfun(@(x) x.optim.symmetry, RLC_tmd_h));
        % Fix ordering
        for i = 1:nS
            if ~isempty(RLC_org_s(i).cross_talk)
                RLC_org(nH+i).cross_talk(2) = RLC_org_s(i).cross_talk(2)+nH;
            end
            if ~isempty(RLC_tmd_s(i).cross_talk)
                RLC_tmd(nH+i).cross_talk(2) = RLC_tmd_s(i).cross_talk(2)+nH;
            end
            RLC_org(nH+i).number            = RLC_org_s(i).number+nH;
            RLC_tmd(nH+i).number            = RLC_tmd_s(i).number+nH;
            RLC_org(nH+i).excitation.entity = RLC_org_s(i).excitation.entity + max_entity_org;
            RLC_tmd(nH+i).excitation.entity = RLC_tmd_s(i).excitation.entity + max_entity_tmd;
            if RLC_tmd(nH+i).optim.boolean == 1
                RLC_tmd(nH+i).optim.symmetry = RLC_tmd_s(i).optim.symmetry + max_symmetry;
            end
        end
        RLC_tmd_w = [];
        RLC_org_w = [];
    % Wire - Shield
    elseif ~isempty(MREDM.WIE.RLC) && isempty(MREDM.SIE.RLC) && ~isempty(MREDM.SIE.shield_RLC)
        RLC_tmd_w = MREDM.WIE.RLC;
        RLC_tmd_h = MREDM.SIE.shield_RLC;
        wire_file = fullfile('./data/coils/wire_files',MREDM.inputs.wire);
        RLC_org_w = geo_scoil_lumped_elements(wire_file,0);
        shield_file = fullfile('./data/coils/shield_files',MREDM.inputs.shield);
        RLC_org_h = geo_scoil_lumped_elements(shield_file,0);
        nH        = RLC_org_h(end).number;
        nW        = RLC_org_w(end).number;
        RLC_org   = repmat(struct('number',[], 'type',[], 'load',[], 'value',[], 'optim',[], 'cross_talk', [], 'excitation', []), nH+nW,1);
        RLC_tmd   = repmat(struct('number',[], 'type',[], 'load',[], 'value',[], 'optim',[], 'cross_talk', [], 'excitation', []), nH+nW,1);
        % Concatenate
        RLC_org(1:nH) = RLC_org_h(:);
        RLC_tmd(1:nH) = RLC_tmd_h(:);
        RLC_org(nH+1:end) = RLC_org_w(:);
        RLC_tmd(nH+1:end) = RLC_tmd_w(:);
        separator_counter = 0;
        for i = 1:length(RLC_tmd_h)
            if RLC_tmd_h(i).optim.boolean == 1
                separator_counter = separator_counter + length(convertCharsToStrings(RLC_tmd_h(i).load));
            end
        end
        max_entity_org = max(arrayfun(@(x) x.excitation.entity, RLC_org_h));
        max_entity_tmd = max(arrayfun(@(x) x.excitation.entity, RLC_tmd_h));
        max_symmetry = max(arrayfun(@(x) x.optim.symmetry, RLC_tmd_h));
        % Fix ordering
        for i = 1:nW
            if ~isempty(RLC_org_w(i).cross_talk)
                RLC_org(nH+i).cross_talk(2) = RLC_org_w(i).cross_talk(2)+nH;
            end
            if ~isempty(RLC_tmd_w(i).cross_talk)
                RLC_tmd(nH+i).cross_talk(2) = RLC_tmd_w(i).cross_talk(2)+nH;
            end
            RLC_org(nH+i).number            = RLC_org_w(i).number+nH;
            RLC_tmd(nH+i).number            = RLC_tmd_w(i).number+nH;
            RLC_org(nH+i).excitation.entity = RLC_org_w(i).excitation.entity + max_entity_org;
            RLC_tmd(nH+i).excitation.entity = RLC_tmd_w(i).excitation.entity + max_entity_tmd;
            if RLC_tmd(nH+i).optim.boolean == 1
                RLC_tmd(nH+i).optim.symmetry = RLC_tmd_w(i).optim.symmetry + max_symmetry;
            end
        end
        RLC_tmd_s = [];
        RLC_org_s = [];
    % Wire - Coil - Shield
    elseif ~isempty(MREDM.WIE.RLC) && ~isempty(MREDM.SIE.RLC) && ~isempty(MREDM.SIE.shield_RLC)
        RLC_tmd_s = MREDM.SIE.RLC;
        RLC_tmd_w = MREDM.WIE.RLC;
        RLC_tmd_h = MREDM.SIE.shield_RLC;
        coil_file = fullfile('./data/coils/coil_files',MREDM.inputs.coil);
        RLC_org_s = geo_scoil_lumped_elements(coil_file,0);
        wire_file = fullfile('./data/coils/wire_files',MREDM.inputs.wire);
        RLC_org_w = geo_scoil_lumped_elements(wire_file,0);
        shield_file = fullfile('./data/coils/shield_files',MREDM.inputs.shield);
        RLC_org_h = geo_scoil_lumped_elements(shield_file,0);
        nS        = RLC_org_s(end).number;
        nH        = RLC_org_h(end).number;
        nW        = RLC_org_w(end).number;
        RLC_org   = repmat(struct('number',[], 'type',[], 'load',[], 'value',[], 'optim',[], 'cross_talk', [], 'excitation', []), nH+nW,1);
        RLC_tmd   = repmat(struct('number',[], 'type',[], 'load',[], 'value',[], 'optim',[], 'cross_talk', [], 'excitation', []), nH+nW,1);
        % Concatenate
        RLC_org(1:nH) = RLC_org_h(:);
        RLC_tmd(1:nH) = RLC_tmd_h(:);
        RLC_org(nH+1:end) = RLC_org_w(:);
        RLC_tmd(nH+1:end) = RLC_tmd_w(:);
        RLC_org(nH+nW+1:end) = RLC_org_s(:);
        RLC_tmd(nH+nW+1:end) = RLC_tmd_s(:);
        separator_counter.first = 0;
        separator_counter.second = 0;
        for i = 1:length(RLC_tmd_h)
            if RLC_tmd_h(i).optim.boolean == 1
                separator_counter.first = separator_counter.first + length(convertCharsToStrings(RLC_tmd_h(i).load));
            end
        end
        for i = 1:length(RLC_tmd_w)
            if RLC_tmd_w(i).optim.boolean == 1
                separator_counter.second = separator_counter.second + length(convertCharsToStrings(RLC_tmd_w(i).load));
            end
        end
        max_entity_org = max(arrayfun(@(x) x.excitation.entity, RLC_org_h));
        max_entity_tmd = max(arrayfun(@(x) x.excitation.entity, RLC_tmd_h));
        max_symmetry   = max(arrayfun(@(x) x.optim.symmetry, RLC_tmd_h));
        % Fix ordering
        for i = 1:nW
            if ~isempty(RLC_org_w(i).cross_talk)
                RLC_org(nH+i).cross_talk(2) = RLC_org_w(i).cross_talk(2)+nH;
            end
            if ~isempty(RLC_tmd_w(i).cross_talk)
                RLC_tmd(nH+i).cross_talk(2) = RLC_tmd_w(i).cross_talk(2)+nH;
            end
            RLC_org(nH+i).number            = RLC_org_w(i).number+nH;
            RLC_tmd(nH+i).number            = RLC_tmd_w(i).number+nH;
            RLC_org(nH+i).excitation.entity = RLC_org_w(i).excitation.entity + max_entity_org;
            RLC_tmd(nH+i).excitation.entity = RLC_tmd_w(i).excitation.entity + max_entity_tmd;
            if RLC_tmd(nH+i).optim.boolean == 1
                RLC_tmd(nH+i).optim.symmetry = RLC_tmd_w(i).optim.symmetry + max_symmetry;
            end
        end
        max_entity_org2 = max(arrayfun(@(x) x.excitation.entity, RLC_org_w));
        max_entity_tmd2 = max(arrayfun(@(x) x.excitation.entity, RLC_tmd_w));
        max_symmetry2   = max(arrayfun(@(x) x.optim.symmetry, RLC_tmd_w));
        % Fix ordering
        for i = 1:nS
            if ~isempty(RLC_org_s(i).cross_talk)
                RLC_org(nH+nW+i).cross_talk(2) = RLC_org_s(i).cross_talk(2)+nH+nW;
            end
            if ~isempty(RLC_tmd_s(i).cross_talk)
                RLC_tmd(nH+nW+i).cross_talk(2) = RLC_tmd_s(i).cross_talk(2)+nH+nW;
            end
            RLC_org(nH+nW+i).number            = RLC_org_s(i).number+nH+nW;
            RLC_tmd(nH+nW+i).number            = RLC_tmd_s(i).number+nH+nW;
            RLC_org(nH+nW+i).excitation.entity = RLC_org_s(i).excitation.entity + max_entity_org + max_entity_org2;
            RLC_tmd(nH+nW+i).excitation.entity = RLC_tmd_s(i).excitation.entity + max_entity_tmd + max_entity_tmd2;
            if RLC_tmd(nH+nW+i).optim.boolean == 1
                RLC_tmd(nH+nW+i).optim.symmetry = RLC_tmd_s(i).optim.symmetry + max_symmetry + max_symmetry2;
            end
        end

    end
end
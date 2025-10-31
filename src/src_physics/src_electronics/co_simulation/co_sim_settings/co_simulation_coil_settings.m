function[RLC_org,RLC_org_w,RLC_org_s,RLC_org_h,wire_file,coil_file] = co_simulation_coil_settings(MREDM)  
    
    wire_file = [];
    coil_file = [];

    % Coil
    if isempty(MREDM.WIE.RLC) && ~isempty(MREDM.SIE.RLC) && isempty(MREDM.SIE.shield_RLC)
        coil_file = fullfile('./data/coils/coil_files',MREDM.inputs.coil);
        RLC_org   = geo_scoil_lumped_elements(coil_file,0);
        RLC_org_s = RLC_org;
        RLC_org_w = [];
        RLC_org_h = [];
    % Wire
    elseif ~isempty(MREDM.WIE.RLC) && isempty(MREDM.SIE.RLC) && isempty(MREDM.SIE.shield_RLC)
        wire_file = fullfile('./data/coils/wire_files',MREDM.inputs.wire);
        RLC_org   = geo_scoil_lumped_elements(wire_file,0);
        RLC_org_w = RLC_org;
        RLC_org_s = [];
        RLC_org_h = [];
    % Shield
    elseif isempty(MREDM.WIE.RLC) && isempty(MREDM.SIE.RLC) && ~isempty(MREDM.SIE.shield_RLC)
        shield_file = fullfile('./data/coils/shield_files',MREDM.inputs.shield);
        RLC_org   = geo_scoil_lumped_elements(shield_file,0);
        RLC_org_h = RLC_org;
        RLC_org_s = [];
        RLC_org_w = [];
    % Wire - Coil
    elseif ~isempty(MREDM.WIE.RLC) && ~isempty(MREDM.SIE.RLC) && isempty(MREDM.SIE.shield_RLC)
        coil_file = fullfile('./data/coils/coil_files',MREDM.inputs.coil);
        RLC_org_s = geo_scoil_lumped_elements(coil_file,0);
        wire_file = fullfile('./data/coils/wire_files',MREDM.inputs.wire);
        RLC_org_w = geo_scoil_lumped_elements(wire_file,0);
        RLC_org_h = [];
        nW        = RLC_org_w(end).number;
        nS        = RLC_org_s(end).number;
        RLC_org   = repmat(struct('number',[], 'type',[], 'load',[], 'value',[], 'optim',[], 'cross_talk', [], 'excitation', []), nW+nS,1);
        % Concatenate
        RLC_org(1:nW) = RLC_org_w(:);
        RLC_org(nW+1:end) = RLC_org_s(:);
        max_entity_org = max(arrayfun(@(x) x.excitation.entity, RLC_org_w));
        % Fix ordering
        for i = 1:nS
            if ~isempty(RLC_org_s(i).cross_talk)
                RLC_org(nW+i).cross_talk(2) = RLC_org_s(i).cross_talk(2)+nW;
            end
            RLC_org(nW+i).number            = RLC_org_s(i).number+nW;
            RLC_org(nW+i).excitation.entity = RLC_org_s(i).excitation.entity + max_entity_org;
        end
    % Coil - Shield
    elseif isempty(MREDM.WIE.RLC) && ~isempty(MREDM.SIE.RLC) && ~isempty(MREDM.SIE.shield_RLC)
        coil_file = fullfile('./data/coils/coil_files',MREDM.inputs.coil);
        RLC_org_s = geo_scoil_lumped_elements(coil_file,0);
        shield_file = fullfile('./data/coils/shield_files',MREDM.inputs.shield);
        RLC_org_h = geo_scoil_lumped_elements(shield_file,0);
        RLC_org_w = [];
        nH        = RLC_org_h(end).number;
        nS        = RLC_org_s(end).number;
        RLC_org   = repmat(struct('number',[], 'type',[], 'load',[], 'value',[], 'optim',[], 'cross_talk', [], 'excitation', []), nH+nS,1);
        % Concatenate
        RLC_org(1:nH) = RLC_org_h(:);
        RLC_org(nH+1:end) = RLC_org_s(:);
        max_entity_org = max(arrayfun(@(x) x.excitation.entity, RLC_org_h));
        % Fix ordering
        for i = 1:nS
            if ~isempty(RLC_org_s(i).cross_talk)
                RLC_org(nH+i).cross_talk(2) = RLC_org_s(i).cross_talk(2)+nH;
            end
            RLC_org(nH+i).number            = RLC_org_s(i).number+nH;
            RLC_org(nH+i).excitation.entity = RLC_org_s(i).excitation.entity + max_entity_org;
        end
    % Wire - Shield
    elseif ~isempty(MREDM.WIE.RLC) && isempty(MREDM.SIE.RLC) && ~isempty(MREDM.SIE.shield_RLC)
        wire_file = fullfile('./data/coils/wire_files',MREDM.inputs.wire);
        RLC_org_w = geo_scoil_lumped_elements(wire_file,0);
        shield_file = fullfile('./data/coils/shield_files',MREDM.inputs.shield);
        RLC_org_h = geo_scoil_lumped_elements(shield_file,0);
        RLC_org_s = [];
        nH        = RLC_org_h(end).number;
        nW        = RLC_org_w(end).number;
        RLC_org   = repmat(struct('number',[], 'type',[], 'load',[], 'value',[], 'optim',[], 'cross_talk', [], 'excitation', []), nW+nH,1);
        % Concatenate
        RLC_org(1:nH) = RLC_org_h(:);
        RLC_org(nH+1:end) = RLC_org_w(:);
        max_entity_org = max(arrayfun(@(x) x.excitation.entity, RLC_org_h));
        % Fix ordering
        for i = 1:nW
            if ~isempty(RLC_org_w(i).cross_talk)
                RLC_org(nH+i).cross_talk(2) = RLC_org_w(i).cross_talk(2)+nH;
            end
            RLC_org(nH+i).number            = RLC_org_w(i).number+nH;
            RLC_org(nH+i).excitation.entity = RLC_org_w(i).excitation.entity + max_entity_org;
        end
    % Wire - Coil - Shield
    elseif ~isempty(MREDM.WIE.RLC) && ~isempty(MREDM.SIE.RLC) && ~isempty(MREDM.SIE.shield_RLC)
        wire_file = fullfile('./data/coils/wire_files',MREDM.inputs.wire);
        RLC_org_w = geo_scoil_lumped_elements(wire_file,0);
        coil_file = fullfile('./data/coils/coil_files',MREDM.inputs.coil);
        RLC_org_s = geo_scoil_lumped_elements(coil_file,0);
        shield_file = fullfile('./data/coils/shield_files',MREDM.inputs.shield);
        RLC_org_h = geo_scoil_lumped_elements(shield_file,0);
        nH        = RLC_org_h(end).number;
        nW        = RLC_org_w(end).number;
        nS        = RLC_org_s(end).number;
        RLC_org   = repmat(struct('number',[], 'type',[], 'load',[], 'value',[], 'optim',[], 'cross_talk', [], 'excitation', []), nW+nS+nH,1);
        % Concatenate
        RLC_org(1:nH) = RLC_org_h(:);
        RLC_org(nH+1:nH+nW) = RLC_org_w(:);
        RLC_org(nH+nW+1:end) = RLC_org_s(:);
        max_entity_org = max(arrayfun(@(x) x.excitation.entity, RLC_org_h));
        % Fix ordering
        for i = 1:nW
            if ~isempty(RLC_org_w(i).cross_talk)
                RLC_org(nH+i).cross_talk(2) = RLC_org_w(i).cross_talk(2)+nH;
            end
            RLC_org(nH+i).number            = RLC_org_w(i).number+nH;
            RLC_org(nH+i).excitation.entity = RLC_org_w(i).excitation.entity + max_entity_org;
        end
        max_entity_org2 = max(arrayfun(@(x) x.excitation.entity, RLC_org_w));
        % Fix ordering
        for i = 1:nS
            if ~isempty(RLC_org_s(i).cross_talk)
                RLC_org(nH+nW+i).cross_talk(2) = RLC_org_s(i).cross_talk(2)+nH+nW;
            end
            RLC_org(nH+nW+i).number            = RLC_org_s(i).number+nH+nW;
            RLC_org(nH+nW+i).excitation.entity = RLC_org_s(i).excitation.entity + max_entity_org+max_entity_org2;
        end
    end
end
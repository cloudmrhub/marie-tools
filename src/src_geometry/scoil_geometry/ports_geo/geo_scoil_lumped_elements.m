function LC = geo_scoil_lumped_elements(filename, tmd)
% GEO_SCOIL_LUMPED_ELEMENTS_JSON
% Reads a JSON-based lumped element definition file for MARIE.
% Recreates the same struct structure and parsing logic as the .txt version.
%
% INPUTS:
%   filename : JSON file path
%   tmd      : tuning/matching/decoupling flag (0/1)
%
% OUTPUT:
%   LC : struct array containing lumped elements
%
% MARIE 3.0
% Author: Ilias I. Giannakopoulos (NYU Grossman School of Medicine, 2025)
% --------------------------------------------------

    % Change file extension
    [filepath, name, ~] = fileparts(filename);
    jsonFile = fullfile(filepath, [name, '.json']);

    % Read and decode JSON
    fid = fopen(jsonFile, 'r');
    raw = fread(fid, inf, 'uint8=>char')';
    fclose(fid);
    data = jsondecode(raw);

    elems = data.coil_configuration.elements;
    N = numel(elems);

    LC = repmat(struct('number', [], 'type', [], 'load', [], 'value', [], ...
                       'optim', [], 'cross_talk', [], 'excitation', []), N, 1);

    for i = 1:N
        e = elems(i);

        % ---------------------------
        % Basic fields
        % ---------------------------
        LC(i).number = e.number;
        LC(i).type   = e.type;
        LC(i).load   = e.load;
        LC(i).optim  = e.optim;
        LC(i).cross_talk = e.cross_talk;
        LC(i).excitation = e.excitation;

        % ---------------------------
        % Handle multi-part matching networks
        % ---------------------------
        if strcmp(e.type, 'port')
            % Split load like "capacitorParallel_inductorSeries_capacitorParallel"
            LC(i).load = strsplit(e.load, '_');

            % Value can be an array of capacitances/inductances
            if iscell(e.value)
                LC(i).value = cellfun(@str2double, e.value);
            else
                LC(i).value = e.value(:)';  % Ensure row vector
            end
        else
            LC(i).load = e.load;
            LC(i).value = e.value;
        end

        % ---------------------------
        % Optimization data
        % ---------------------------
        if strcmp(e.type, 'port')
            % Port has vector min/max arrays and symmetry
            LC(i).optim.boolean = e.optim.boolean;
            LC(i).optim.minim   = e.optim.minim(:)'; % ensure row
            LC(i).optim.maxim   = e.optim.maxim(:)'; % ensure row
            LC(i).optim.symmetry = e.optim.symmetry;
        else
            LC(i).optim.boolean = e.optim.boolean;
            LC(i).optim.minim   = e.optim.minim;
            LC(i).optim.maxim   = e.optim.maxim;
            LC(i).optim.symmetry = e.optim.symmetry;

            % Convert eligible elements to ports in TMD mode
            if tmd && LC(i).optim.boolean
                LC(i).type = 'port';
            end
        end

        % ---------------------------
        % Cross-talk (for mutual inductors)
        % ---------------------------
        if any(strcmp(LC(i).load, 'mutual_inductor')) && isstruct(e.cross_talk)
            if isfield(e.cross_talk, 'coupled_port')
                LC(i).cross_talk = [];
                LC(i).cross_talk(1) = e.cross_talk.coupled_port;
                LC(i).cross_talk(2) = e.cross_talk.coupled_value;
            end
        elseif isempty(e.cross_talk)
            LC(i).cross_talk = [];
        end

        % ---------------------------
        % Excitation
        % ---------------------------
        LC(i).excitation.entity = e.excitation.entity;
        LC(i).excitation.TxRx   = e.excitation.TxRx;

        % Voltage logic
        if strcmp(LC(i).type, 'port')
            LC(i).excitation.voltage = 1;
        elseif strcmp(LC(i).type, 'element') && strcmp(LC(i).load, 'port')
            LC(i).excitation.voltage = 1;
        else
            LC(i).excitation.voltage = 0;
        end
    end
end

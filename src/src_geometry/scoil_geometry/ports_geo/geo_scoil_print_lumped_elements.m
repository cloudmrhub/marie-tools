function N_ports = geo_scoil_print_lumped_elements(filename, LC, LC_tmd, new_elems)
    
    [filepath, name, ~] = fileparts(filename);
    filename = fullfile(filepath, [name, '.json']);

    fid = fopen(filename, 'r');
    if fid == -1, error('Cannot open input JSON: %s', filename); end
    raw = fread(fid, inf, 'uint8=>char')';
    fclose(fid);
    orig = jsondecode(raw);

    if ~isfield(orig, 'coil_configuration') || ~isfield(orig.coil_configuration, 'elements')
        error('Unexpected JSON layout: missing coil_configuration.elements');
    end

    c1 = 0;
    N_ports = 0;

    for i = 1:numel(LC)
        LC(i).optim = LC_tmd(i).optim;

        if LC_tmd(i).optim.boolean == 1
            if strcmp(LC(i).type, 'port')
                for ii = 1:numel(LC(i).value)
                    c1 = c1 + 1;
                    LC(i).value(ii) = new_elems(c1);
                end
            elseif strcmp(LC(i).type, 'element')
                c1 = c1 + 1;
                LC(i).value = new_elems(c1);
            end
        end

        if strcmp(LC(i).type, 'port')
            N_ports = N_ports + 1;
        end
    end

    elems_out = LC;
    for i = 1:numel(elems_out)
        if strcmp(elems_out(i).type, 'port')
            if iscell(elems_out(i).load)
                elems_out(i).load = strjoin(elems_out(i).load, '_'); 
            elseif isstring(elems_out(i).load)
                elems_out(i).load = char(elems_out(i).load);
            end
            if isfield(elems_out(i).optim,'minim'), elems_out(i).optim.minim = elems_out(i).optim.minim(:)'; end
            if isfield(elems_out(i).optim,'maxim'), elems_out(i).optim.maxim = elems_out(i).optim.maxim(:)'; end
            elems_out(i).value = elems_out(i).value(:)'; % vector for ports
        else
            if iscell(elems_out(i).load)
                elems_out(i).load = strjoin(elems_out(i).load, '_');
            end
        end

        if isfield(elems_out(i),'excitation') && isfield(elems_out(i).excitation,'voltage')
            elems_out(i).excitation = rmfield(elems_out(i).excitation,'voltage');
        end

        if isempty(elems_out(i).cross_talk)
            elems_out(i).cross_talk = struct();
        end
    end

    orig.coil_configuration.elements = elems_out;

    fid = fopen(filename, 'w');
    if fid == -1, error('Cannot open output JSON for writing: %s', filename); end
    json_str = jsonencode(orig, 'PrettyPrint', true);
    fwrite(fid, json_str, 'char');
    fclose(fid);
end

function [] = to_CloudMR(input_file,solution,author_name)
    
    %% Set Paths
    input_file    = fullfile('./data/inputs/', input_file);
    solution_file = fullfile('./data/solutions/',solution);
    storer_dir    = strcat(erase(solution,'.mat'));
    folderName    = './data/solutions/CloudMR_hub/data';
    FigureDirName = './data/solutions/CloudMR_hub/other';
    fileName      = './data/solutions/CloudMR_hub/info.json';
    hub_dir       = './data/solutions/CloudMR_hub/';
    CoilFolder    = './data/solutions/CloudMR_hub/data/coil/';

    reset_folder(folderName);
    reset_folder(FigureDirName);
    reset_folder(CoilFolder);

    %% Read Inputs
    fid  = fopen(input_file, 'r');
    raw  = fread(fid, inf, 'uint8=>char')';
    fclose(fid);
    inputs_json = jsondecode(raw);
    
    B0            = inputs_json.B0;                   
    NUCLEUS       = inputs_json.Nucleus;            
    RHBM          = inputs_json.BodyFile;                    
    COIL          = inputs_json.CoilFile;                    
    WIRE          = inputs_json.WireFile;                    
    SHIELD        = inputs_json.ShieldFile;                  
    BASIS_SUPPORT = inputs_json.SurfaceBasisSupportFile;                   
    BASIS_FILE    = inputs_json.BasisFile;            
    if double(inputs_json.Basis_Functions_VIE) == 0
        VIE_BASIS = 'PWC';
    elseif double(inputs_json.Basis_Functions_VIE) == 1
        VIE_BASIS = 'PWL';
    end                  
    if double(inputs_json.TMD) == 1
        TMD = 'YES';
    elseif double(inputs_json.TMD) == 0
        TMD = 'NO';
    end       
    [~, cmdout] = system('grep "model name" /proc/cpuinfo | uniq');
    cpuName = strtrim(cmdout); 
    
    %% Settings
    S = load(strcat('./data/bodies/',RHBM));
    res_x  = squeeze(abs(S.RHBM.r(2,1,1,1)-S.RHBM.r(1,1,1,1)));
    res    = [res_x,res_x,res_x];
    origin = squeeze(S.RHBM.r(1,1,1,:));

    %% DATA
    phi = h5read_or_zero(solution_file,'/solution/e_field/Rx/phi') + h5read_or_zero(solution_file,'/solution/e_field/Rx/phi_coil') + h5read_or_zero(solution_file,'/solution/e_field/Rx/phi_lump');

    b1m = h5read_or_zero(solution_file,'/solution/h_field/Rx/b1m');

    b1p = h5read_or_zero(solution_file,'/solution/h_field/Tx/b1p');

    Tx_ports = abs(isscalar(b1p)-1)*size(b1p,4);
    Rx_ports = abs(isscalar(b1m)-1)*size(b1m,4);

    emc = em_constants(B0,NUCLEUS);

    epsilon_r = field_or_zero(S.RHBM,'epsilon_r');
    sigma_e   = field_or_zero(S.RHBM,'sigma_e');
    RhoH      = field_or_zero(S.RHBM,'rho');
    mu_r      = field_or_zero(S.RHBM,'mu_r');
    c         = field_or_zero(S.RHBM,'c');
    k         = field_or_zero(S.RHBM,'k');
    w         = field_or_zero(S.RHBM,'w');
    Q_m       = field_or_zero(S.RHBM,'Q_m');
    T1        = field_or_zero(S.RHBM,'T1');
    T2        = field_or_zero(S.RHBM,'T2');
    RhoM      = field_or_zero(S.RHBM,'mass_rho');
    Dw        = field_or_zero(S.RHBM,'dw');
    T2star    = field_or_zero(S.RHBM,'T2s');

    %% OUTPUT FILE
    MARIE_OUT.headers.platform.username = getenv('USER');
    MARIE_OUT.headers.platform.matlab_version = version;
    MARIE_OUT.headers.platform.os = system_dependent('getos');
    MARIE_OUT.headers.platform.computer_name = computer;
    MARIE_OUT.headers.platform.gpu = gpuDevice().Name;
    MARIE_OUT.headers.platform.cpu = cpuName;
    MARIE_OUT.headers.version = '3.0.1';
    MARIE_OUT.headers.author = author_name;

    MARIE_OUT.headers.Inputs.b0 = B0;
    MARIE_OUT.headers.Inputs.nucleus = NUCLEUS;
    MARIE_OUT.headers.Inputs.freq = emc.freq;
    MARIE_OUT.headers.Inputs.basis_vie = VIE_BASIS;
    MARIE_OUT.headers.Inputs.basis_sie = 'RWG';
    MARIE_OUT.headers.Inputs.basis_wie = 'Triangle';
    MARIE_OUT.headers.Inputs.phantom = 'Duke_2mm';
    MARIE_OUT.headers.Inputs.Number_of_tissues = length(unique(sigma_e))-1;
    MARIE_OUT.headers.Inputs.resolution = res;
    MARIE_OUT.headers.Inputs.coil = COIL;
    MARIE_OUT.headers.Inputs.wire = WIRE;
    MARIE_OUT.headers.Inputs.shield = SHIELD;
    MARIE_OUT.headers.Inputs.basis_support = BASIS_SUPPORT;
    MARIE_OUT.headers.Inputs.basis_file = BASIS_FILE;  
    MARIE_OUT.headers.Inputs.Number_of_Rx_channels = Rx_ports;
    MARIE_OUT.headers.Inputs.Number_of_Tx_channels = Tx_ports;
    MARIE_OUT.headers.Inputs.co_simulation = TMD;
    MARIE_OUT.headers.Inputs.EM_Simulator = 'MARIE_3.0_WSVIE_version';

    %% JSON File
    maps = {
    T1,        'T1',                    't1',        'data/t1.nii.gz',        1000;
    T2,        'T2',                    't2',        'data/t2.nii.gz',        1001;
    T2star,    'T2_star',               't2star',    'data/t2star.nii.gz',    1003;
    RhoH,      'Proton_Density',        'rhoh',      'data/rhoh.nii.gz',      1004;
    RhoM,      'Mass_Density',          'rhom',      'data/rhom.nii.gz',      1005;
    Dw,        'Chemical_Shift',        'dw',        'data/dw.nii.gz',        1006;
    epsilon_r, 'Relative_Permittivity', 'epsilonr',  'data/epsilon_r.nii.gz', 1007;
    sigma_e,   'Conductivity',          'sigmae',    'data/sigma_e.nii.gz',   1008;
    mu_r,      'Relative_Permeability', 'mu_r',      'data/mur.nii.gz',       1009;
    c,         'Heat_Capacity',         'c',         'data/c.nii.gz',         1010;
    k,         'Thermal_Conductivity',  'k',         'data/k.nii.gz',         1011;
    w,         'Perfusion',             'w',         'data/w.nii.gz',         1012;
    Q_m,       'Heat Generation Rate',  'Q_m',       'data/q.nii.gz',         1013;};


    assets = {
    'Coil mesh',     'coilmsh',   'coil_files',   COIL,                 'coil',  true;
    'Wire mesh',     'wiremsh',   'wire_files',   WIRE,                 'wire',  true;
    'Shield mesh',   'shieldmsh', 'shield_files', SHIELD,               'shield',true;
    'Basis support', 'basismsh',  'basis_files',  BASIS_SUPPORT,        'basis', false;};


    g = struct([]);
    for i = 1:size(maps,1)
        data = maps{i,1};
        if isequal(data,0) || isempty(data)
            continue
        end
        g(i).name  = maps{i,2};
        g(i).description = maps{i,3};
        g(i).filename = maps{i,4};
        g(i).id = maps{i,5};
        g(i).type ='accessory';
        g(i).matPixelType = class(data);
        g(i).pixelType = 'real';
        g(i).dim = 3;
        g(i).numberOfComponentsPerPixel = 1;
        save_nii( ...
            make_nii(squeeze(data),res,origin), ...
            fullfile(hub_dir, g(i).filename) ...
        );
    end

    g_assets = empty_entry();
    g_assets = g_assets([]);   
    
    for i = 1:size(assets,1)
        name  = assets{i,1};
        desc  = assets{i,2};
        fold  = assets{i,3};
        file  = assets{i,4};
        id    = assets{i,5};
        has_json = assets{i,6};
        mesh_src = fullfile('data','coils',fold, file);
        mesh_dst = fullfile(CoilFolder,replace(file,'/','_'));
    
        copy_if_exists(mesh_src, mesh_dst);
    
        g_assets(end+1) = build_asset(name, desc, fullfile('data','coils',fold,file), id);
        if has_json
            ports = replace(file,'.msh','.json');  
            ports_src = fullfile('data','coils',fold, ports);
            ports_dst = fullfile(CoilFolder,replace(ports,'/','_'));
            copy_if_exists(ports_src, ports_dst);
            g_assets(end+1) = build_asset( ...
                [name ' ports'], ...
                [desc 'ports'], ...
                fullfile('data/coils/coil_files',ports), ...
                [id '_ports']);
        end
    end

    g_phi = write_output(phi, 'Noise Covariance Matrix', 'noisecovariancematrix', 'data/psi.nii.gz', 1, 2, hub_dir, res, origin, false);

    g_b1m = write_channels(b1m, 'Magnetic Receive Field', 'b1m', 99, hub_dir, res, origin);

    g_b1p = write_channels(b1p, 'Magnetic Transmit Field', 'b1p', 299, hub_dir, res, origin);

    %% Visualize
    coil = []; wire = []; shield = [];
    if ~isempty(COIL)
        [coil.node,~,~,coil.elem,~] = Mesh_Parse(strcat('./data/coils/coil_files/',COIL));
    end
    if ~isempty(WIRE)
        [~,~,~,wire.F_point,wire.S_point,~] = Mesh_Wire(strcat('./data/coils/wire_files/',WIRE));
    end
    if ~isempty(SHIELD)
        [shield.node,~,~,shield.elem,~] = Mesh_Parse(strcat('./data/coils/shield_files/',SHIELD));
    end
    visualize_body_wire_coil_shield(S.RHBM,wire,coil,shield)
    FigureName = fullfile(FigureDirName, 'Geometry.jpg');
    saveas(gcf, FigureName);
    
    %% STORE DATA    
    MARIE_OUT.data = [g_phi,g_b1m,g_b1p,g_assets,g];
    
    jsonData = jsonencode(MARIE_OUT);
    fileID = fopen(fileName, 'w');
    fprintf(fileID, '%s', jsonData);
    fclose(fileID);
    
    cd(hub_dir)
    zipFilename = strcat(storer_dir,'.zip');
    zip(zipFilename, {'./data', 'info.json','./other'});
    
    cd('../../../')    
    delete(fileName);
    rmdir(folderName, 's');
    rmdir(FigureDirName, 's');
end

function data = h5read_or_zero(file, path)
    try
        h5info(file, path);   
        data = h5read(file, path);
        data = data.real+1i*data.imag;
    catch
        data = 0;            
    end
end
function out = field_or_zero(S, field)
    if isfield(S, field) && ~isempty(S.(field))
        out = rot90(S.(field),2);
    else
        out = 0;
    end
end
function g = build_asset(name, desc, filename, id)
    g.name = name;
    g.description = desc;
    g.filename = filename;
    g.id = id;
    g.type ='accessory';
    g.matPixelType = NaN;
    g.pixelType = NaN;
    g.dim = NaN;
    g.numberOfComponentsPerPixel = NaN;
end
function copy_if_exists(src, dst)
    if isfile(src)
        copyfile(src, dst);
    end
end
function g = write_output(data, name, desc, filename, id, dim, hub_dir, res, origin, rotate)
    if isempty(data) || (isscalar(data) && isequal(data,0))
        g = [];    
        return
    end
    g = empty_entry(); 
    if nargin < 10
        rotate = false;
    end
    if rotate
        data = rot90(data,2);
    end
    g.name = name;
    g.description = desc;
    g.filename = filename;
    g.id = id;
    g.type ='output';
    g.matPixelType = class(data);
    g.pixelType = 'complex';
    g.dim = dim;
    g.numberOfComponentsPerPixel = 1;
    nii = make_nii(squeeze(data),res,origin);
    save_nii(nii, fullfile(hub_dir,filename));
end
function g = write_channels(vol4D, prefix, desc, id_offset, hub_dir, res, origin)
    if isempty(vol4D) || (isscalar(vol4D) && vol4D==0)
        g = empty_entry();
        g = g([]);   
        return
    end
    nCh = size(vol4D,4);
    tmp = cell(1,nCh);
    for i = 1:nCh
        N = sprintf('%03d',i);
        tmp{i} = write_output( ...
            vol4D(:,:,:,i), ...
            sprintf('%s %d',prefix,i), ...
            desc, ...
            sprintf('data/%s_%s.nii.gz',desc,N), ...
            id_offset + i, ...
            3, ...
            hub_dir, res, origin, true);
    end
    g = [tmp{:}]; 
end

function g = empty_entry()
    g = struct( ...
        'name', [], ...
        'description', [], ...
        'filename', [], ...
        'id', [], ...
        'type', [], ...
        'matPixelType', [], ...
        'pixelType', [], ...
        'dim', [], ...
        'numberOfComponentsPerPixel', [] );
end
function reset_folder(folder)
    if ~isfolder(folder)
        mkdir(folder);
        return
    end
    files = dir(fullfile(folder,'*'));
    files = files(~[files.isdir]); 
    for k = 1:numel(files)
        delete(fullfile(folder,files(k).name));
    end
    dirs = dir(folder);
    dirs = dirs([dirs.isdir]);
    dirs = dirs(~ismember({dirs.name},{'.','..'}));
    for k = 1:numel(dirs)
        rmdir(fullfile(folder,dirs(k).name),'s');
    end
end



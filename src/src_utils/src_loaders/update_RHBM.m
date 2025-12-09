function [MREDM] = update_RHBM(MREDM)

    % Create Full Domain
    temp = load(MREDM.inputs.rhbm);    
    RHBM = temp.RHBM;
    RHBM.name = 'Master_Domain';
    RHBM.epsilon_r = 10*ones(size(RHBM.epsilon_r));
    RHBM.sigma_e = 0.1*ones(size(RHBM.sigma_e));
    RHBM.idxS = find(RHBM.sigma_e>0);
    RHBM.rho = ones(size(RHBM.rho));

    % Update File's Name
    idx = find(MREDM.inputs.rhbm == '/', 1, 'last');
    MREDM.inputs.rhbm = [MREDM.inputs.rhbm(1:idx) RHBM.name '.mat'];

    % Save Master File
    save(MREDM.inputs.rhbm, 'RHBM');

end
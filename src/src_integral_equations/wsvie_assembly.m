function [MREDM] = wsvie_assembly(MREDM,geo_flag) 

    fprintf('\tAssembling WIE-Coil Operators ...\n');
    MREDM = wie_assembly(MREDM);
    fprintf('\tAssembling SIE-Coil Operators ...\n');
    MREDM = sie_assembly(MREDM);
    fprintf('\tAssembling SIE-Shield Operators ...\n');
    MREDM = shield_sie_assembly(MREDM);
    fprintf('\tAssembling VIE-Body Operators ...\n');
    MREDM = vie_assembly(MREDM,geo_flag);
    fprintf('\tAssembling WSVIE-Coupling Operators ...\n'); 
    MREDM = wsvie_coupling_assembly(MREDM);

end
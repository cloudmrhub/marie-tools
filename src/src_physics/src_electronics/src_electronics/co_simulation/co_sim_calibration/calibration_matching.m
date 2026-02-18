function [M_trans] = calibration_matching(Ss,Sm,z0)

    road = 'schur';
    
    n = length(Ss);
    I = eye(n);
    
    Sa22  = Ss';
    if strcmp(road,'eig')
        [V,D] = eig(I-Ss'*Ss);
        Sa21  = V*sqrt(D)/V;
    elseif strcmp(road,'schur')
        [Q,U] = schur(I-Ss'*Ss);
        Sa21 = Q*sqrt(U)*Q';
    end
    Sa12  = Sa21.';
    Sa11  = -(Sa12')\Ss*Sa21;
    Sa    = [Sa11 Sa12; Sa21 Sa22];

    Sb11  = Sm;
    if strcmp(road,'eig')
        [V,D] = eig(I-Sm'*Sm);
        Sb21  = V*sqrt(D)/V;
    elseif strcmp(road,'schur')
        [Q,U] = schur(I-Sm'*Sm);
        Sb21 = Q*sqrt(U)*Q';
    end
    Sb12  = Sb21.';
    Sb22  = ((Sb12'*Sb11)/Sb21)';
    Sb    = [Sb11 Sb12; Sb21 Sb22];

    abcd1 = np_s2abcd(Sb,z0);
    abcd2 = np_s2abcd(Sa,z0);

    Sc    = np_abcd2s(abcd1*abcd2,z0); 
    Sc21  = Sc(n+1:end,1:n);
    Sc22  = Sc(n+1:end,n+1:end);   

    M_trans=(I-Sc22*Ss)\Sc21;

    % Sc11  = Sc(1:n,1:n);
    % Sc12  = Sc(1:n,n+1:end);
    % S_tran = Sc11+Sc12*Ss*((I-Sc22*Ss)\Sc21);

end
